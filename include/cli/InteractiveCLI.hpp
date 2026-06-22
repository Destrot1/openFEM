#pragma once
#include "fem/Model.hpp"
#include "geometry/StepImporter.hpp"
#include "mesh/Mesher.hpp"
#include <string>

class InteractiveCLI {
public:
    void run();

private:
    FEMModel     m_model;
    StepImporter m_importer;
    Mesher       m_mesher;

    void cmdLoad();
    void cmdAssign();
    void cmdMesh();
    void cmdBoundary();
    void cmdLoads();
    void cmdSolve();
    void cmdResults();
    void cmdSave();
    void cmdLoadJson();
    void cmdSnapshot();
    void cmdRestore();
    void cmdRefresh();

    // Re-exports scene.json after any command that changes the model,
    // so the browser viewer (which polls for changes) picks it up on
    // its own -- there's no user-facing 'web' command anymore. Silent
    // no-op if there's no geometry yet.
    void autoExportWeb();

    void        printBanner();
    void        printStatus();
    void        printMenu();
    double      promptDouble(const std::string& msg, double def = 0.0);
    std::string promptString(const std::string& msg,
                             const std::string& def = "");
    bool        promptConfirm(const std::string& msg);
};