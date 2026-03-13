#include "model.h"
#include "reader.h"
#include "solver.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>

static void print_usage(const char *progname) {
    fprintf(stderr,
            "Usage: %s <scp_file> [options]\n"
            "  --verbosity N       Verbosity level (default 2)\n"
            "  --max-nodes N       Maximum BnB nodes (default 100000)\n"
            "  --time-limit S      Time limit in seconds (default 0 = unlimited)\n"
            "  --cut-frequency F   Cut frequency 0.0-1.0 (default 0.0, 0=disable)\n"
            "  --balas-frequency F Balas frequency 0.0-1.0 (default 0.6, 0=disable)\n"
            "  --cut-rounds N      Mid-BnB cut rounds per event (default 3)\n"
            "  --branch STRATEGY   Branch strategy: reliability, most_fractional, highest_cost (default reliability)\n"
            "  --reliability-eta N Reliability parameter (default 4)\n"
            "  --reliability-sb N  Max strong-branch probes per node (default 8)\n"
            "  --show-solution     Print selected columns\n"
            "  --preprocess RULES  Preprocess rules (default \"single,two\")\n",
            progname);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string input_file;
    scpsol::SolverConfig config;
    bool found_file = false;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--verbosity") == 0 && i + 1 < argc) {
            config.verbosity = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--max-nodes") == 0 && i + 1 < argc) {
            config.max_nodes = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--time-limit") == 0 && i + 1 < argc) {
            config.time_limit_seconds = atof(argv[++i]);
        } else if (strcmp(argv[i], "--cut-frequency") == 0 && i + 1 < argc) {
            config.mid_bnb_cut_frequency = atof(argv[++i]);
        } else if (strcmp(argv[i], "--balas-frequency") == 0 && i + 1 < argc) {
            config.aggressive_balas_frequency = atof(argv[++i]);
        } else if (strcmp(argv[i], "--cut-rounds") == 0 && i + 1 < argc) {
            config.mid_bnb_cut_rounds = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--show-solution") == 0) {
            config.show_solution = true;
        } else if (strcmp(argv[i], "--branch") == 0 && i + 1 < argc) {
            config.branch_strategy = argv[++i];
        } else if (strcmp(argv[i], "--reliability-eta") == 0 && i + 1 < argc) {
            config.reliability_eta = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--reliability-sb") == 0 && i + 1 < argc) {
            config.reliability_max_sb = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--preprocess") == 0 && i + 1 < argc) {
            config.preprocess_rules = argv[++i];
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        } else if (argv[i][0] != '-' && !found_file) {
            input_file = argv[i];
            found_file = true;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    if (!found_file) {
        fprintf(stderr, "Error: no input file specified\n");
        print_usage(argv[0]);
        return 1;
    }

    // Derive enable flags from frequency parameters
    config.cuts_enabled = (config.mid_bnb_cut_frequency > 0.0);
    config.balas_enabled = (config.aggressive_balas_frequency > 0.0);

    try {
        scpsol::ScpInstance instance = scpsol::read_scp_file(input_file);
        if (config.verbosity >= 1) {
            fprintf(stderr, "Instance: %d rows, %d cols, %d nnz\n",
                    instance.nrows, instance.ncols,
                    static_cast<int>(instance.csr_values.size()));
        }

        scpsol::SolverResult result = scpsol::solve(instance, config);

        // Output
        printf("Status:    %s\n", result.status.c_str());
        printf("Primal:    %.12g\n", result.primal_obj);
        printf("Dual:      %.12g\n", result.dual_obj);
        if (std::isfinite(result.mip_gap))
            printf("MIP gap:   %.6f%%\n", result.mip_gap * 100.0);
        else
            printf("MIP gap:   inf\n");
        printf("Nodes:     %d\n", result.nodes_processed);
        printf("LP solves: %d\n", result.lp_solves);
        printf("Wall time: %.3f s\n", result.wall_time);

        if (config.show_solution && !result.solution.empty()) {
            printf("Solution:\n");
            for (int j = 0; j < static_cast<int>(result.solution.size()); ++j) {
                if (result.solution[static_cast<size_t>(j)] > 0.5)
                    printf("  x[%d] = 1\n", j);
            }
        }

        return (result.status == "Optimal" || result.status == "TimeLimit" ||
                result.status == "NodeLimit") ? 0 : 1;
    } catch (const std::exception &e) {
        fprintf(stderr, "Error: %s\n", e.what());
        return 1;
    }
}
