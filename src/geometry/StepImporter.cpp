#include "geometry/StepImporter.hpp"
#include <fmt/core.h>

#include <STEPControl_Reader.hxx>
#include <BRepTools.hxx>

// g_shape is a file-scope global -- the importer is effectively a
// singleton and cannot hold more than one geometry at a time
static TopoDS_Shape g_shape;

bool StepImporter::load(const std::string& filepath) {
    fmt::print("[StepImporter] Loading: {}\n", filepath);

    STEPControl_Reader reader;
    if (reader.ReadFile(filepath.c_str()) != IFSelect_RetDone) {
        fmt::print("[StepImporter] ERROR: cannot read the file.\n");
        return false;
    }
    reader.TransferRoots();
    g_shape = reader.OneShape();
    return true;
}

bool StepImporter::exportToBRep(const std::string& path) const {
    BRepTools::Write(g_shape, path.c_str());
    fmt::print("[StepImporter] BREP written to: {}\n", path);
    return true;
}
