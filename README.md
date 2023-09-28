Overview of project stages:

```mermaid
graph TD
classDef Title stroke-width:0, color:blue, font-weight:bold, font-size: 19px;

0[Functions of <br> ConsensusOPLS-DA method]
    subgraph Number1[Step 0 : Translate MATLAB function into R]
    0 --> 1
    1(Step0_1_matrix2saisir) 
    0 --> 2
    2(RVConsensusOPLSPerm)

    1 --- |include in matrix2saisir| 3
    3(add code)

    2 --> 4(RVConsensusOPLS)
    2 --> 5(RV_modified)

    4 --> 6
    6[koplsScale]
    4 --> 7
    7[koplsKernel]
    4 --> 8
    8[koplsModel]
    4 -->|use| 5
    4 --> 9
    9[ConsensusOPLSCV]
    4 --> 10
    10[DQ2]

    9 --- |see graph below| 9
    end
    class Number1 Title;
```

To make it easier to read, we have reproduced the previous graph using the 
ConsensusOPLSCV function only.

```mermaid
graph TB
classDef Title stroke-width:0, color:blue, font-weight:bold, font-size: 19px;

0[Functions of <br> ConsensusOPLS-DA method]
    subgraph Number1[Step 0 : Translate MATLAB function into R <br> from the ConsensusOPLSCV]
    0 --> 1
    1(ConsensusOPLSCV) 
	
	1 --> 2
    2(koplsDummy)
    1 --> 3
    3(koplsReDummy)
    1 --> 4
    4(koplsCrossValSet)
    1 --> 5
    5(koplsScale)
    1 --> 6
    6(koplsScaleApply)
    1 --> 7
    7(koplsCenterKT*)
    1 ---|use koplsModel| 1
    1 --> 8
    8(koplsPredict)
    1 --> 9
    9(koplsRescale)
    1 --> 10
    10(koplsMaxClassify)
    1 --> 11
    11(koplsBasicClassify)
    1 --> 12
    12(koplsSensSpec)
    1 --> 13
    13(koplsConfusionMatrix)

    end
    class Number1 Title;
```

This structure has been reproduced for the git tree: one branch per function 
and per code file. Next, the functions were tested on example datasets and 
verified. Finally, all the branches will be merged to finalize the method.

Once the method had been translated, verified and validated, it was tested 
on a real dataset.

(Coming soon)

The next step is to introduce variable selection to this ConsensusOPLS-DA 
method.

(Coming soon)