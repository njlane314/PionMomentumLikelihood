#include <iostream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include <torch/torch.h>

struct PionData {
    double sim_px;
    double sim_py;
    double sim_pz;

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

    double sim_px, sim_py, sim_pz;
    std::vector<double> *sp_x = nullptr, *sp_y = nullptr, *sp_z = nullptr, *sp_e = nullptr;

    tree->SetBranchAddress("px", &sim_px);
    tree->SetBranchAddress("py", &sim_py);
    tree->SetBranchAddress("pz", &sim_pz);
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

        PionData pionData;
        pionData.sim_px = sim_px;
        pionData.sim_py = sim_py;
        pionData.sim_pz = sim_pz;

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
    std::vector<double> momenta;
    for (const auto& entry : data) {
        auto x = vectorToTensor(entry.sp_x);
        auto y = vectorToTensor(entry.sp_y);
        auto z = vectorToTensor(entry.sp_z);
        auto e = vectorToTensor(entry.sp_e);

        auto traj = torch::stack({x, y, z, e}, 1);  
        traj = padToMaxLength(traj, max_length); 
        tensors.push_back(traj);
        momenta.push_back(entry.sim_px);
        momenta.push_back(entry.sim_py);
        momenta.push_back(entry.sim_pz);
    }

    if (tensors.empty()) {
        throw std::runtime_error("No valid data found in the ROOT file.");
    }

    auto states = torch::stack(tensors);
    auto labels = torch::tensor(momenta, torch::kFloat32).view({-1, 3});

    return {states, labels};
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
        fc = register_module("fc", torch::nn::Linear(d_model, 3));  
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

        auto dataset = CustomDataset(x, y).map(torch::data::transforms::Stack<>());
        auto dataloader = torch::data::make_data_loader(dataset, 64);

        Transformer model(d_model, 4, 2);
        torch::optim::Adam optimizer(model->parameters(), torch::optim::AdamOptions(0.001));

        for (size_t epoch = 0; epoch < 10; ++epoch) {
            torch::Tensor loss;
            for (auto& batch : *dataloader) {
                auto inputs = batch.data;
                auto targets = batch.target;

                optimizer.zero_grad();
                auto output = model->forward(inputs);
                loss = torch::mse_loss(output, targets);
                loss.backward();
                optimizer.step();
            }
            std::cout << "Epoch [" << epoch << "], Loss: " << loss.item<double>() << std::endl;
        }

        torch::NoGradGuard no_grad;
        model->eval();
        torch::Tensor outputs, targets;
        for (auto& batch : *dataloader) {
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

    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}