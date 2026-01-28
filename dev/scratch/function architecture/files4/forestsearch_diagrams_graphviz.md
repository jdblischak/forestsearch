# ForestSearch Architecture - DiagrammeR/Graphviz Version
# 
# This file contains Graphviz DOT code that can be rendered using:
# 1. R: DiagrammeR::grViz() function
# 2. Online: https://dreampuf.github.io/GraphvizOnline/
# 3. Command line: dot -Tpng file.dot -o file.png
#
# Each diagram is in a separate code block below.

## How to use in R:

```r
# Install DiagrammeR if needed
# install.packages("DiagrammeR")

library(DiagrammeR)

# Copy the DOT code between the triple backticks and paste into grViz()
grViz("
  digraph {
    ...
  }
")
```

---

## Diagram 1: High-Level Pipeline

```dot
digraph ForestSearchPipeline {
  rankdir=LR;
  node [shape=box, style="rounded,filled", fontname="Arial"];
  
  // Nodes
  input [label="Input\nData", fillcolor="#E3F2FD"];
  varsel [label="Variable\nSelection", fillcolor="#E8F4FD"];
  dataprep [label="Data\nPreparation", fillcolor="#FFF3E0"];
  search [label="Subgroup\nSearch", fillcolor="#FCE4EC"];
  consist [label="Consistency\nEvaluation", fillcolor="#E8F5E9"];
  output [label="Output", fillcolor="#F3E5F5"];
  
  // Edges
  input -> varsel -> dataprep -> search -> consist -> output;
}
```

---

## Diagram 2: Main Function Flow

```dot
digraph ForestSearchMain {
  rankdir=TB;
  node [shape=box, style="rounded,filled", fontname="Arial"];
  
  // Main function
  fs [label="forestsearch()", shape=ellipse, fillcolor="#2E86AB", fontcolor="white", style="filled"];
  
  // Steps
  p1 [label="1. Validate Inputs", fillcolor="#ECEFF1"];
  p2 [label="2. GRF Variable Selection\ngrf.subg.harm.survival()", fillcolor="#E8F4FD"];
  p3 [label="3. LASSO Selection\nlasso_selection()", fillcolor="#E8F4FD"];
  p4 [label="4. Data Preparation\nget_FSdata()", fillcolor="#FFF3E0"];
  p5 [label="5. Subgroup Search\nsubgroup.search()", fillcolor="#FCE4EC"];
  p6 [label="6. Consistency Evaluation\nsubgroup.consistency()", fillcolor="#E8F5E9"];
  p7 [label="7. Return Results", fillcolor="#F3E5F5"];
  
  // Flow
  fs -> p1 -> p2 -> p3 -> p4 -> p5 -> p6 -> p7;
}
```

---

## Diagram 3: GRF Variable Selection

```dot
digraph GRFSelection {
  rankdir=TB;
  node [shape=box, style="rounded,filled", fontname="Arial", fillcolor="#E8F4FD"];
  
  input [label="Input Data\n+ Confounders", fillcolor="#E3F2FD"];
  fit [label="fit_causal_forest()"];
  vi [label="Compute Variable\nImportance"];
  tree [label="Fit Policy Trees\n(depth 1, 2, 3)"];
  select [label="select_best_subgroup()"];
  extract [label="extract_all_tree_cuts()"];
  output [label="GRF Cuts\n(e.g., 'age > 50')", fillcolor="#E8F5E9"];
  
  input -> fit -> vi -> tree -> select -> extract -> output;
}
```

---

## Diagram 4: LASSO Selection

```dot
digraph LASSOSelection {
  rankdir=TB;
  node [shape=box, style="rounded,filled", fontname="Arial", fillcolor="#E8F4FD"];
  
  input [label="Confounders", fillcolor="#E3F2FD"];
  cv [label="cv.glmnet()\nCox-LASSO"];
  lambda [label="Select lambda.min"];
  coef [label="Extract Coefficients"];
  
  // Decision
  decision [label="Coef ≠ 0?", shape=diamond, fillcolor="#FFFDE7"];
  selected [label="SELECTED", fillcolor="#E8F5E9"];
  omitted [label="OMITTED", fillcolor="#FFEBEE"];
  
  input -> cv -> lambda -> coef -> decision;
  decision -> selected [label="Yes"];
  decision -> omitted [label="No"];
}
```

---

## Diagram 5: Data Preparation

```dot
digraph DataPrep {
  rankdir=TB;
  node [shape=box, style="rounded,filled", fontname="Arial", fillcolor="#FFF3E0"];
  
  input [label="Raw Data", fillcolor="#E3F2FD"];
  classify [label="Classify Variables\nis.continuous()"];
  
  // Branch
  cont [label="Continuous\nCreate Cut Indicators"];
  cat [label="Categorical\ndummy()"];
  
  lasso [label="Apply LASSO Filter\nfilter_by_lassokeep()"];
  grf [label="Apply GRF Cuts"];
  force [label="Handle Forced Cuts\nget_conf_force()"];
  output [label="Binary Indicator\nMatrix Z", fillcolor="#E8F5E9"];
  
  input -> classify;
  classify -> cont;
  classify -> cat;
  cont -> lasso;
  cat -> lasso;
  lasso -> grf -> force -> output;
}
```

---

## Diagram 6: Subgroup Search Algorithm

```dot
digraph SubgroupSearch {
  rankdir=TB;
  node [shape=box, style="rounded,filled", fontname="Arial"];
  
  gen [label="Generate All\nCombinations", fillcolor="#E3F2FD"];
  loop [label="Loop: Each\nCombination", fillcolor="#FFF3E0"];
  
  // Checks
  prev [label="Prevalence ≥ minp?", shape=diamond, fillcolor="#FFFDE7"];
  size [label="Size ≥ n.min?", shape=diamond, fillcolor="#FFFDE7"];
  events [label="Events OK?\nd0 ≥ d0.min\nd1 ≥ d1.min", shape=diamond, fillcolor="#FFFDE7"];
  
  cox [label="Fit Cox Model\nCompute HR", fillcolor="#E8F4FD"];
  hrcheck [label="HR > threshold?", shape=diamond, fillcolor="#FFFDE7"];
  store [label="Store Candidate", fillcolor="#E8F5E9"];
  output [label="Return Sorted\nCandidates", fillcolor="#E8F5E9"];
  
  gen -> loop -> prev;
  prev -> size [label="Yes"];
  prev -> loop [label="No"];
  size -> events [label="Yes"];
  size -> loop [label="No"];
  events -> cox [label="Yes"];
  events -> loop [label="No"];
  cox -> hrcheck;
  hrcheck -> store [label="Yes"];
  hrcheck -> loop [label="No"];
  store -> loop;
  loop -> output [label="Done"];
}
```

---

## Diagram 7: Fixed Consistency Algorithm

```dot
digraph FixedConsistency {
  rankdir=TB;
  node [shape=box, style="rounded,filled", fontname="Arial"];
  
  input [label="Candidate\nSubgroup", fillcolor="#E3F2FD"];
  loop [label="Repeat n.splits times", fillcolor="#FFF3E0"];
  split [label="Random 50/50 Split", fillcolor="#E8F4FD"];
  hr1 [label="HR in Split 1"];
  hr2 [label="HR in Split 2"];
  check [label="Both HR >\nthreshold?", shape=diamond, fillcolor="#FFFDE7"];
  yes [label="Consistent\nCount++", fillcolor="#E8F5E9"];
  no [label="Inconsistent", fillcolor="#FFEBEE"];
  result [label="Pcons =\nConsistent / Total", fillcolor="#E8F5E9"];
  
  input -> loop -> split;
  split -> hr1;
  split -> hr2;
  hr1 -> check;
  hr2 -> check;
  check -> yes [label="Yes"];
  check -> no [label="No"];
  yes -> loop;
  no -> loop;
  loop -> result [label="Done"];
}
```

---

## Diagram 8: Two-Stage Consistency Algorithm

```dot
digraph TwoStageConsistency {
  rankdir=TB;
  node [shape=box, style="rounded,filled", fontname="Arial"];
  
  input [label="Candidate", fillcolor="#E3F2FD"];
  stage1 [label="Stage 1: Quick Screen\n(30 splits)", fillcolor="#E8F4FD"];
  check1 [label="Pass threshold?", shape=diamond, fillcolor="#FFFDE7"];
  reject1 [label="REJECT", fillcolor="#FFEBEE"];
  
  stage2 [label="Stage 2: Batched\nEvaluation", fillcolor="#FFF3E0"];
  batch [label="Run 20 splits"];
  wilson [label="Compute Wilson CI"];
  checklow [label="CI lower >\nthreshold?", shape=diamond, fillcolor="#FFFDE7"];
  accept [label="ACCEPT", fillcolor="#E8F5E9"];
  checkup [label="CI upper <\nthreshold?", shape=diamond, fillcolor="#FFFDE7"];
  reject2 [label="REJECT", fillcolor="#FFEBEE"];
  
  input -> stage1 -> check1;
  check1 -> reject1 [label="No"];
  check1 -> stage2 [label="Yes"];
  stage2 -> batch -> wilson -> checklow;
  checklow -> accept [label="Yes"];
  checklow -> checkup [label="No"];
  checkup -> reject2 [label="Yes"];
  checkup -> stage2 [label="No\n(continue)"];
}
```

---

## Diagram 9: Bootstrap Bias Correction

```dot
digraph Bootstrap {
  rankdir=TB;
  node [shape=box, style="rounded,filled", fontname="Arial"];
  
  input [label="Original Result\n(H_obs)", fillcolor="#E3F2FD"];
  gen [label="Generate B\nBootstrap Samples", fillcolor="#E8F4FD"];
  loop [label="For each bootstrap b", fillcolor="#FFF3E0"];
  draw [label="Draw bootstrap\nsample"];
  run [label="Run forestsearch()"];
  compute [label="Compute:\nH*, H**, H*obs"];
  agg [label="Aggregate Results", fillcolor="#E8F5E9"];
  bias [label="Compute Bias\nCorrection", fillcolor="#E8F5E9"];
  output [label="Adjusted HR", fillcolor="#E8F5E9"];
  
  input -> gen -> loop -> draw -> run -> compute;
  compute -> loop [label="next b"];
  loop -> agg [label="Done"];
  agg -> bias -> output;
}
```

---

## Diagram 10: K-Fold Cross-Validation

```dot
digraph KFold {
  rankdir=TB;
  node [shape=box, style="rounded,filled", fontname="Arial"];
  
  input [label="Full Dataset", fillcolor="#E3F2FD"];
  split [label="Split into K Folds", fillcolor="#E8F4FD"];
  loop [label="For fold k = 1 to K", fillcolor="#FFF3E0"];
  train [label="Train: K-1 folds"];
  test [label="Test: fold k"];
  fs [label="forestsearch()\non training"];
  eval [label="Evaluate on test"];
  agg [label="Aggregate Results", fillcolor="#E8F5E9"];
  
  input -> split -> loop;
  loop -> train;
  loop -> test;
  train -> fs -> eval;
  test -> eval;
  eval -> loop [label="next k"];
  loop -> agg [label="Done"];
}
```

---

## Diagram 11: Package Dependencies

```dot
digraph Dependencies {
  rankdir=TB;
  node [shape=box, style="rounded,filled", fontname="Arial"];
  
  fs [label="forestsearch", fillcolor="#2E86AB", fontcolor="white"];
  
  // Statistical
  subgraph cluster_stat {
    label="Core Statistical";
    style=filled;
    fillcolor="#E8F4FD";
    survival [label="survival"];
    grf [label="grf"];
    glmnet [label="glmnet"];
    policytree [label="policytree"];
  }
  
  // Data
  subgraph cluster_data {
    label="Data Handling";
    style=filled;
    fillcolor="#FFF3E0";
    dt [label="data.table"];
    stringr [label="stringr"];
  }
  
  // Parallel
  subgraph cluster_par {
    label="Parallel Computing";
    style=filled;
    fillcolor="#FCE4EC";
    future [label="future"];
    doFuture [label="doFuture"];
    foreach [label="foreach"];
    callr [label="callr"];
  }
  
  // Visualization
  subgraph cluster_viz {
    label="Visualization";
    style=filled;
    fillcolor="#E8F5E9";
    ggplot2 [label="ggplot2"];
    gt [label="gt"];
    forestploter [label="forestploter"];
  }
  
  fs -> survival;
  fs -> grf;
  fs -> glmnet;
  fs -> policytree;
  fs -> dt;
  fs -> stringr;
  fs -> future;
  fs -> doFuture;
  fs -> foreach;
  fs -> callr;
  fs -> ggplot2;
  fs -> gt;
  fs -> forestploter;
}
```

---

## How to Render in R

```r
library(DiagrammeR)

# Example: Render the high-level pipeline
grViz('
digraph ForestSearchPipeline {
  rankdir=LR;
  node [shape=box, style="rounded,filled", fontname="Arial"];
  
  input [label="Input\\nData", fillcolor="#E3F2FD"];
  varsel [label="Variable\\nSelection", fillcolor="#E8F4FD"];
  dataprep [label="Data\\nPreparation", fillcolor="#FFF3E0"];
  search [label="Subgroup\\nSearch", fillcolor="#FCE4EC"];
  consist [label="Consistency\\nEvaluation", fillcolor="#E8F5E9"];
  output [label="Output", fillcolor="#F3E5F5"];
  
  input -> varsel -> dataprep -> search -> consist -> output;
}
')

# Save to file
diagram <- grViz('...')
export_svg(diagram) %>% charToRaw() %>% rsvg::rsvg_png("diagram.png")
```

---

## Online Rendering

Copy any DOT code block (without the ```dot markers) to:
- https://dreampuf.github.io/GraphvizOnline/
- https://edotor.net/
- https://viz-js.com/
