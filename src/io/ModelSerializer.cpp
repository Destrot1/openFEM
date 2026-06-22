#include "io/ModelSerializer.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <fmt/core.h>

using json = nlohmann::json;

bool ModelSerializer::save(const FEMModel& model, const std::string& path) {
    json j;
    j["name"]     = model.name;
    j["stepFile"] = model.stepFile;

    // Materials
    for (const auto& m : model.materials) {
        j["materials"].push_back({
            {"id",            m.id},
            {"name",          m.name},
            {"E",             m.E},
            {"nu",            m.nu},
            {"rho",           m.rho},
            {"yieldStrength", m.yieldStrength}
        });
    }

    // Beam sections
    for (const auto& s : model.sections) {
        j["sections"].push_back({
            {"id",   s.id},
            {"name", s.name},
            {"area", s.area},
            {"Iy",   s.Iy},
            {"Iz",   s.Iz},
            {"J",    s.J}
        });
    }

    // Assignments
    for (const auto& a : model.assignments) {
        j["assignments"].push_back({
            {"dim",        a.dim},
            {"entityId",   a.entityId},
            {"femType",    (int)a.femType},
            {"materialId", a.materialId},
            {"thickness",  a.thickness},
            {"sectionId",  a.sectionId},
            {"meshSize",   a.meshSize}
        });
    }

    // Boundary conditions
    for (const auto& bc : model.bcs) {
        j["bcs"].push_back({
            {"dim",      bc.dim},
            {"entityId", bc.entityId},
            {"type",     (int)bc.type},
            {"locked",   std::vector<bool>(bc.locked, bc.locked+6)}
        });
    }

    // Loads
    for (const auto& l : model.loads) {
        j["loads"].push_back({
            {"dim",      l.dim},
            {"entityId", l.entityId},
            {"type",     (int)l.type},
            {"values",   std::vector<double>(l.values, l.values+3)}
        });
    }

    std::ofstream f(path);
    if (!f) {
        fmt::print("[Serializer] ERROR: cannot write {}\n", path);
        return false;
    }
    f << j.dump(2);
    fmt::print("[Serializer] Model saved to: {}\n", path);
    return true;
}

bool ModelSerializer::load(FEMModel& model, const std::string& path) {
    std::ifstream f(path);
    if (!f) {
        fmt::print("[Serializer] ERROR: file not found {}\n", path);
        return false;
    }

    json j;
    f >> j;

    model.name     = j.value("name",     "");
    model.stepFile = j.value("stepFile", "");

    model.materials.clear();
    for (const auto& m : j.value("materials", json::array())) {
        model.materials.push_back({
            m["id"], m["name"], m["E"], m["nu"], m["rho"],
            m.value("yieldStrength", 0.0)
        });
    }

    model.sections.clear();
    for (const auto& s : j.value("sections", json::array())) {
        BeamSection sec;
        sec.id   = s["id"];
        sec.name = s["name"];
        sec.area = s["area"];
        sec.Iy   = s["Iy"];
        sec.Iz   = s["Iz"];
        sec.J    = s["J"];
        model.sections.push_back(sec);
    }

    model.assignments.clear();
    for (const auto& a : j.value("assignments", json::array())) {
        EntityAssignment ea;
        ea.dim        = a["dim"];
        ea.entityId   = a["entityId"];
        ea.femType    = (ElementType)(int)a["femType"];
        ea.materialId = a["materialId"];
        ea.thickness  = a["thickness"];
        ea.sectionId  = a["sectionId"];
        ea.meshSize   = a["meshSize"];
        model.assignments.push_back(ea);
    }

    model.bcs.clear();
    for (const auto& bc : j.value("bcs", json::array())) {
        BoundaryCondition c;
        c.dim      = bc["dim"];
        c.entityId = bc["entityId"];
        c.type     = (BCType)(int)bc["type"];
        auto locks = bc["locked"].get<std::vector<bool>>();
        for (int i = 0; i < 6 && i < (int)locks.size(); ++i)
            c.locked[i] = locks[i];
        model.bcs.push_back(c);
    }

    model.loads.clear();
    for (const auto& l : j.value("loads", json::array())) {
        Load ld;
        ld.dim      = l["dim"];
        ld.entityId = l["entityId"];
        ld.type     = (LoadType)(int)l["type"];
        auto vals   = l["values"].get<std::vector<double>>();
        for (int i = 0; i < 3 && i < (int)vals.size(); ++i)
            ld.values[i] = vals[i];
        model.loads.push_back(ld);
    }

    model.geometryLoaded  = !model.stepFile.empty();
    model.physicsAssigned = !model.assignments.empty();

    fmt::print("[Serializer] Model loaded from: {}\n", path);
    return true;
}
