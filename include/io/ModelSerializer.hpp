#pragma once
#include "fem/Model.hpp"
#include <string>

class ModelSerializer {
public:
    static bool save(const FEMModel& model, const std::string& path);
    static bool load(FEMModel& model,       const std::string& path);
};