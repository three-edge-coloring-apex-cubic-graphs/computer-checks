#include "cartwheel.hpp"
#include "configuration.hpp"
#include "configuration_homomorphism.hpp"
#include "pseudo_configuration.hpp"
#include "pseudo_embedding.hpp"
#include "rule.hpp"
#include <boost/program_options.hpp>
#include <filesystem>
#include <fstream>
#include <spdlog/spdlog.h>

using namespace boost::program_options;
namespace fs = std::filesystem;

bool existArgs(const variables_map &vm, const vector<string> &args) {
    for (const auto &arg : args) {
        if (!vm.count(arg)) {
            spdlog::warn(fmt::format("Specify {}.", arg));
            return false;
        }
    }
    return true;
}

optional<ConfForm> valid_extension(const string &extension) {
    if (extension == "outer-extension") {
        return ConfForm::OuterExtension;
    } else if (extension == "free-completion") {
        return ConfForm::FreeCompletion;
    } else {
        spdlog::warn(fmt::format("Invalid extension: {}.", extension));
        return std::nullopt;
    }
}

int main(const int ac, const char *const *const av) {
    using namespace boost::program_options;
    options_description description("Options");
    // clang-format off
    description.add_options()
        ("combine_rules", "Combine rules")
        ("enum_wheels", "Enumerate wheels")
        ("verify_no_bad_cartwheels", "Verify that no bad cartwheels exist")
        ("all_hom_images", "Enumerate all homomorphic images of a configuration with no obstructing cycle")
        ("confdir,C", value<string>(), "A directory containing configuration files")
        ("ruledir,R", value<string>(), "A directory containing rule files")
        ("combined_ruledir,S", value<string>(), "A directory containing combined rule files")
        ("rule_auxiliary_dir,A", value<string>(), "A directory containing rule auxiliary files")
        ("cartwheeldir,W", value<string>(), "A directory containing cartwheel files")
        ("wheel,w", value<string>(), "A wheel file")
        ("degree,d", value<int>(), "Degree of the center vertex of wheels")
        ("max_deg2_vertices", value<int>(), "Maximum number of degree-2 vertices in a multi-boundary island")
        ("extension", value<string>(), "Extension kind (outer-extension, free-completion)")
        ("outdir,o", value<string>(), "Output directory")
        ("help,H", "Display options")
        ("verbosity,v", value<int>()->default_value(0), "1 for debug, 2 for trace");
    // clang-format on

    variables_map vm;
    store(parse_command_line(ac, av, description), vm);
    notify(vm);

    if (vm.count("help")) {
        description.print(std::cout);
        return 0;
    }
    if (vm.count("verbosity")) {
        auto v = vm["verbosity"].as<int>();
        if (v == 1) {
            spdlog::set_level(spdlog::level::debug);
        }
        if (v == 2) {
            spdlog::set_level(spdlog::level::trace);
        }
    }

    if (vm.count("combine_rules")) {
        if (!existArgs(vm, {"ruledir", "confdir", "outdir"})) {
            exit(1);
        }
        string ruledir = vm["ruledir"].as<string>();
        string confdir = vm["confdir"].as<string>();
        string outdir  = vm["outdir"].as<string>();
        run_combine_rules(confdir, ruledir, outdir);
    }
    if (vm.count("enum_wheels")) {
        if (!existArgs(vm, {"degree", "confdir", "ruledir", "combined_ruledir", "outdir"})) {
            exit(1);
        }
        int degree              = vm["degree"].as<int>();
        string confdir          = vm["confdir"].as<string>();
        string ruledir          = vm["ruledir"].as<string>();
        string combined_ruledir = vm["combined_ruledir"].as<string>();
        string outdir           = vm["outdir"].as<string>();
        run_enum_wheels(degree, confdir, ruledir, combined_ruledir, outdir);
    }
    if (vm.count("verify_no_bad_cartwheels")) {
        if (!existArgs(vm, {"wheel", "confdir", "ruledir", "combined_ruledir", "rule_auxiliary_dir",
                            "outdir"})) {
            exit(1);
        }
        string wheel_file         = vm["wheel"].as<string>();
        string confdir            = vm["confdir"].as<string>();
        string ruledir            = vm["ruledir"].as<string>();
        string combined_ruledir   = vm["combined_ruledir"].as<string>();
        string rule_auxiliary_dir = vm["rule_auxiliary_dir"].as<string>();
        string outdir             = vm["outdir"].as<string>();
        run_verify_no_bad_cartwheels(wheel_file, confdir, ruledir, combined_ruledir,
                                     rule_auxiliary_dir, outdir);
    }
    if (vm.count("all_hom_images")) {
        if (!existArgs(vm, {"confdir", "outdir", "max_deg2_vertices", "extension"})) {
            exit(1);
        }
        string confdir               = vm["confdir"].as<string>();
        string outdir                = vm["outdir"].as<string>();
        int max_deg2_vertices        = vm["max_deg2_vertices"].as<int>();
        string extension             = vm["extension"].as<string>();
        optional<ConfForm> conf_form = valid_extension(extension);
        if (!conf_form) {
            exit(1);
        }
        all_hom_images(confdir, outdir, max_deg2_vertices, *conf_form);
    }
    return 0;
}