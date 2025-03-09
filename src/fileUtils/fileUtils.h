#pragma once

#include <filesystem>

struct FileUtils
{
    static bool findTheFolder(const std::string &sName, std::filesystem::path &path);
};
