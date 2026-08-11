# Computer checks

This repository contains programs for computing free homomorphic images, computing combined rules, enumerating cartwheels, and recursively computing homomorphic images of configurations.

Each lemma in the paper can be verified by running the following commands.

## Preparation
For preparation, we build the program and make directories.
```bash
cmake -S . -B build
cmake --build build
mkdir -p log/d{7..11} empty combined_rules/non_blocked wheels/d{7..11} wheels/positive
```

Clone the repositories containing configuration files and discharging rule files, and auxiliary rule files.
```bash
git clone git@github.com:three-edge-coloring-apex-cubic-graphs/configurations.git
git clone git@github.com:three-edge-coloring-apex-cubic-graphs/discharging-rules.git
```

## Lemma B.1 (for Lemma 6.5)
```bash
./build/src/main --combine_rules -R discharging-rules/R -C configurations/K -o combined_rules/non_blocked > log/non_blocked.log
```

## Lemma B.2 (for Lemma 6.6)
First, we enumerate wheels by considering all possible degrees of the center and its neighbors in the following command.
```bash
bash enum_possible_bad_wheels.sh
```

For each degree $d \in \{7,8,9,10,11\}$, after enumerating wheels of center degree $d$, we enumerate bad cartwheels from these wheels.
In the first step, the number of resulting wheels is decided, so we have only to execute the following commands.
```bash
bash verify_no_bad_cartwheels_for_all.sh 7 4438
bash verify_no_bad_cartwheels_for_all.sh 8 4939
bash verify_no_bad_cartwheels_for_all.sh 9 2409
bash verify_no_bad_cartwheels_for_all.sh 10 567
bash verify_no_bad_cartwheels_for_all.sh 11 38
```

## Lemma B.3  (for Lemma 7.5)
The following command constructs the set of multi-boudary islands $\mathcal{I}$ from the set of configurations $\mathcal{K}$.
```bash
./build/src/main --all_hom_images -C configurations/K -o I --max_deg2_vertices 3 --extension outer-extension > log/I.log &
```
Check the semi-reducibility of all the resulting multi-boundary islands using [the semi-reducibility checker](https://github.com/three-edge-coloring-apex-cubic-graphs/semi-reducibility-checker).
