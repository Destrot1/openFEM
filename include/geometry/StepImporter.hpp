#pragma once
#include <string>

class StepImporter {
public:
    bool load(const std::string& filepath);
    bool exportToBRep(const std::string& outputPath) const;
};
