#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include <torch/torch.h>

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>

struct PionData {
    double sim_px;
    double sim_py;
    double sim_pz;
    
    double sim_purity; 
    double sim_completeness;

    std::vector<double> sp_x;
    std::vector<double> sp_y;
    std::vector<double> sp_z;
    std::vector<double> sp_e;
};

std::vector<PionData> readRootFile(const std::string& fileName) {
    std::vector<PionData> pionDataList;

    TFile *file = TFile::Open(fileName.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return pionDataList;
    }

    TTree *tree = (TTree*)file->Get("ana/Tree");

    if (!tree) {
        std::cerr << "Error: Trees not found in file!" << std::endl;
        return pionDataList;
    }

    double sim_px, sim_py, sim_pz, sim_purity, sim_completeness;
    std::vector<double> *sp_x = nullptr, *sp_y = nullptr, *sp_z = nullptr, *sp_e = nullptr;

    tree->SetBranchAddress("px", &sim_px);
    tree->SetBranchAddress("py", &sim_py);
    tree->SetBranchAddress("pz", &sim_pz);
    tree->SetBranchAddress("purity", &sim_purity); 
    tree->SetBranchAddress("completeness", &sim_completeness);
    tree->SetBranchAddress("x", &sp_x);
    tree->SetBranchAddress("y", &sp_y);
    tree->SetBranchAddress("z", &sp_z);
    tree->SetBranchAddress("e", &sp_e);

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        if (sp_x->empty() || sp_y->empty() || sp_z->empty() || sp_e->empty()) {
            std::cout << "Skipping entry " << i << " due to empty vectors.\n";
            continue; 
        }

        double vec_x = sp_x->back() - sp_x->front();
        double vec_y = sp_y->back() - sp_y->front();
        double vec_z = sp_z->back() - sp_z->front();

        double dot_product = vec_x * sim_px + vec_y * sim_py + vec_z * sim_pz;

        if (dot_product < 0) {
            std::reverse(sp_x->begin(), sp_x->end());
            std::reverse(sp_y->begin(), sp_y->end());
            std::reverse(sp_z->begin(), sp_z->end());
            std::reverse(sp_e->begin(), sp_e->end());
        }

        PionData pionData;
        pionData.sim_px = sim_px;
        pionData.sim_py = sim_py;
        pionData.sim_pz = sim_pz;
        pionData.sim_purity = sim_purity;
        pionData.sim_completeness = sim_completeness;

        pionData.sp_x = *sp_x;
        pionData.sp_y = *sp_y;
        pionData.sp_z = *sp_z;
        pionData.sp_e = *sp_e;

        pionDataList.push_back(pionData);
    }

    file->Close();
    delete file;

    return pionDataList;
}

torch::Tensor vectorToTensor(const std::vector<double>& vec) {
    return torch::tensor(vec, torch::kFloat32);
}

torch::Tensor padToMaxLength(const torch::Tensor& tensor, int64_t max_length) {
    int64_t seq_len = tensor.size(0);
    if (seq_len < max_length) {
        return torch::constant_pad_nd(tensor, {0, 0, 0, max_length - seq_len}, 0);
    } else {
        return tensor.slice(0, 0, max_length); 
    }
}

std::pair<torch::Tensor, torch::Tensor> pionDataToTensor(const std::vector<PionData>& data, int64_t max_length) {
    std::vector<torch::Tensor> tensors;
    std::vector<double> labels;
    for (const auto& entry : data) {
        auto x = vectorToTensor(entry.sp_x);
        auto y = vectorToTensor(entry.sp_y);
        auto z = vectorToTensor(entry.sp_z);
        auto e = vectorToTensor(entry.sp_e);

        auto traj = torch::stack({x, y, z, e}, 1);  
        traj = padToMaxLength(traj, max_length); 
        tensors.push_back(traj);
        labels.push_back(entry.sim_px);
        labels.push_back(entry.sim_py);
        labels.push_back(entry.sim_pz);
        labels.push_back(entry.sim_purity);
        labels.push_back(entry.sim_completeness);
    }

    if (tensors.empty()) {
        throw std::runtime_error("No valid data found in the ROOT file.");
    }

    auto states = torch::stack(tensors);
    auto labels_tensor = torch::tensor(labels, torch::kFloat32).view({-1, 5}); 

    return {states, labels_tensor};
}

class CustomDataset : public torch::data::datasets::Dataset<CustomDataset> {
private:
    torch::Tensor states_, labels_;

public:
    CustomDataset(const torch::Tensor& states, const torch::Tensor& labels)
        : states_(states), labels_(labels) {}

    torch::data::Example<> get(size_t index) override {
        return {states_[index], labels_[index]};
    }

    torch::optional<size_t> size() const override {
        return states_.size(0);
    }
};

struct TransformerImpl : torch::nn::Module {
    TransformerImpl(int64_t d_model, int64_t nhead, int64_t num_encoder_layers)
        : transformer_layer(torch::nn::TransformerEncoderLayerOptions(d_model, nhead)),
          transformer(torch::nn::TransformerEncoderOptions(transformer_layer, num_encoder_layers)) {
        register_module("transformer", transformer);
        fc = register_module("fc", torch::nn::Linear(d_model, 5));  
    }

    torch::Tensor forward(torch::Tensor src) {
        auto output = transformer(src);
        return fc(output.mean(1));
    }

    torch::nn::TransformerEncoderLayer transformer_layer{nullptr};
    torch::nn::TransformerEncoder transformer{nullptr};
    torch::nn::Linear fc{nullptr};
};
TORCH_MODULE(Transformer);

double computeCorrelation(const torch::Tensor& x, const torch::Tensor& y) {
    double mean_x = x.mean().item<double>();
    double mean_y = y.mean().item<double>();

    auto diff_x = x - mean_x;
    auto diff_y = y - mean_y;

    double numerator = (diff_x * diff_y).sum().item<double>();
    double denominator = std::sqrt((diff_x * diff_x).sum().item<double>() * (diff_y * diff_y).sum().item<double>());

    return numerator / denominator;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root>" << std::endl;
        return 1;
    }

    std::string fileName = argv[1];
    const int64_t max_length = 100;

    try {
        auto [x, y] = pionDataToTensor(readRootFile(fileName), max_length);

        const int64_t d_model = 4;

        x = x.view({x.size(0), max_length, d_model});

        auto dataset_size = x.size(0);
        auto train_size = static_cast<int64_t>(dataset_size * 0.8);
        auto test_size = dataset_size - train_size;

        auto train_indices = torch::randperm(dataset_size).narrow(0, 0, train_size);
        auto test_indices = torch::randperm(dataset_size).narrow(0, train_size, test_size);

        auto train_dataset = CustomDataset(x.index_select(0, train_indices), y.index_select(0, train_indices)).map(torch::data::transforms::Stack<>());
        auto test_dataset = CustomDataset(x.index_select(0, test_indices), y.index_select(0, test_indices)).map(torch::data::transforms::Stack<>());

        auto train_dataloader = torch::data::make_data_loader(train_dataset, 64);
        auto test_dataloader = torch::data::make_data_loader(test_dataset, 64);

        Transformer model(d_model, 4, 2);
        torch::optim::Adam optimizer(model->parameters(), torch::optim::AdamOptions(0.001));

        std::vector<double> epoch_losses;

        for (size_t epoch = 0; epoch < 10; ++epoch) {
            model->train();
            torch::Tensor loss;
            for (auto& batch : *train_dataloader) {
                auto inputs = batch.data;
                auto targets = batch.target;

                optimizer.zero_grad();
                auto output = model->forward(inputs);
                loss = torch::mse_loss(output, targets);
                loss.backward();
                optimizer.step();
            }
            std::cout << "Epoch [" << epoch << "], Loss: " << loss.item<double>() << std::endl;
            epoch_losses.push_back(loss.item<double>());
        }

        torch::NoGradGuard no_grad;
        model->eval();
        torch::Tensor outputs, targets;
        for (auto& batch : *test_dataloader) {
            auto inputs = batch.data;
            auto batch_targets = batch.target;

            auto batch_outputs = model->forward(inputs);
            if (!outputs.defined()) {
                outputs = batch_outputs;
                targets = batch_targets;
            } else {
                outputs = torch::cat({outputs, batch_outputs}, 0);
                targets = torch::cat({targets, batch_targets}, 0);
            }
        }

        auto mse = torch::mse_loss(outputs, targets);
        auto mae = torch::mean(torch::abs(outputs - targets));

        std::cout << "Final MSE: " << mse.item<double>() << std::endl;
        std::cout << "Final MAE: " << mae.item<double>() << std::endl;

        // ROOT Plotting
        TCanvas *c1 = new TCanvas("c1", "Training Loss", 800, 600);
        TH1F *h1 = new TH1F("h1", "Training Loss per Epoch;Epoch;Loss", 10, 0, 10);
        for (size_t i = 0; i < epoch_losses.size(); ++i) {
            h1->SetBinContent(i+1, epoch_losses[i]);
        }
        h1->Draw();

        TCanvas *c2 = new TCanvas("c2", "Predictions vs Targets", 800, 600);
        TGraph *g1 = new TGraph(outputs.size(0));
        for (int i = 0; i < outputs.size(0); ++i) {
            g1->SetPoint(i, targets[i][0].item<double>(), outputs[i][0].item<double>());
        }
        g1->SetTitle("Predictions vs Targets;Targets;Predictions");
        g1->Draw("AP");

        // Plot individual targets
        std::vector<const char*> target_labels = {"px", "py", "pz", "purity", "completeness"};
        std::vector<TCanvas*> canvases;
        std::vector<TGraph*> graphs;

        for (int t = 0; t < 5; ++t) {
            TCanvas *c = new TCanvas(("c" + std::to_string(t + 3)).c_str(), ("Predictions vs Targets " + std::string(target_labels[t])).c_str(), 800, 600);
            TGraph *g = new TGraph(outputs.size(0));
            for (int i = 0; i < outputs.size(0); ++i) {
                g->SetPoint(i, targets[i][t].item<double>(), outputs[i][t].item<double>());
            }
            g->SetTitle((std::string("Predictions vs Targets ") + target_labels[t] + ";Targets;Predictions").c_str());
            g->Draw("AP");
            canvases.push_back(c);
            graphs.push_back(g);
        }

        // Plot histograms of prediction errors
        std::vector<TH1F*> error_histograms;
        for (int t = 0; t < 5; ++t) {
            TH1F *h = new TH1F(("h_error_" + std::to_string(t)).c_str(), (std::string("Prediction Error ") + target_labels[t] + ";Error;Counts").c_str(), 100, -1, 1);
            for (int i = 0; i < outputs.size(0); ++i) {
                h->Fill(outputs[i][t].item<double>() - targets[i][t].item<double>());
            }
            error_histograms.push_back(h);
        }

        TCanvas *c8 = new TCanvas("c8", "Prediction Errors", 1600, 1200);
        c8->Divide(3, 2);
        for (int t = 0; t < 5; ++t) {
            c8->cd(t + 1);
            error_histograms[t]->Draw();
        }

        // Plot correlation matrix
        TCanvas *c9 = new TCanvas("c9", "Correlation Matrix", 800, 600);
        TH2F *h_corr = new TH2F("h_corr", "Correlation Matrix;Targets;Predictions", 5, 0, 5, 5, 0, 5);
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                double correlation = computeCorrelation(outputs.select(1, i), targets.select(1, j));
                h_corr->SetBinContent(i + 1, j + 1, correlation);
            }
        }
        h_corr->Draw("COLZ");

        c1->SaveAs("training_loss.png");
        c2->SaveAs("predictions_vs_targets.png");
        for (int t = 0; t < 5; ++t) {
            canvases[t]->SaveAs((std::string("predictions_vs_targets_") + target_labels[t] + ".png").c_str());
        }
        c8->SaveAs("prediction_errors.png");
        c9->SaveAs("correlation_matrix.png");

    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
