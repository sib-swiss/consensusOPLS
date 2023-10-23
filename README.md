# Overview of project stages:

The aim of this project is to translate the `ConsensusOPLS` method from its 
original MATLAB version (available on the GitLab 
https://gitlab.unige.ch/Julien.Boccard/consensusopls) to an R version.

For this purpose, a mind map of the method's functions was created:

```mermaid
graph TD
classDef Title stroke-width:0, color:blue, font-weight:bold, font-size: 19px;

0[Functions of ConsensusOPLS-DA method]
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

0[Functions of ConsensusOPLS-DA method]
    subgraph Number1[Step 0 : from the ConsensusOPLSCV]
    0 --> 1
    1(ConsensusOPLSCV) 
	
	1 --> 2
    2(koplsDummy)
    1 --> 3
    3(koplsReDummy)
    1 --> 4
    4(koplsScaleApply)
    1 --> 5
    5(koplsCenterKT*)
    1 --> 6
    6(koplsPredict)
    1 --> 7
    7(koplsRescale)
    1 --> 8
    8(koplsMaxClassify)
    1 --> 9
    9(koplsBasicClassify)
    1 --> 10
    10(koplsSensSpec)
    1 --> 11
    11(koplsConfusionMatrix)
    1 --> 12
    12(koplsCrossValSet)
    6 --> 13
    13(koplsRescale)

    1 ---|use koplsModel <br> koplsScale | 1
    12 --> 3
    5 ---|KTeTe, KTeTr, KTrTr|5
    6 --> 5


    end
    class Number1 Title;
```

This structure has been reproduced for the git tree: one branch per function 
and per code file. Next, the functions were tested on example datasets 
(`demo_data` proposed by Julien Boccard in its Matlab version) and verified. 
Finally, all the branches will be merged to finalize the method.

# KOPLS R package

At the meeting on 09/10/2023, it turned out that the KOPLS package had already 
been translated into R. To avoid redundancy or errors, the source code of the 
kopls package was searched for. A version was found on the following GitHub: 
https://github.com/sdechaumet/ramopls/tree/master/inst/package (version of 
05/08/2020). Codes previously translated from Matlab to R were compared with 
those in the GitHub repository. References to the authors of the KOPLS package 
have been added at the start of each function's code.

# Partnership with Switzerland

The `main` branch of this project is shared with members of the Swiss Institute 
of Bioinformatics (SIB) on the following GitHub: 
https://github.com/sib-swiss/consensusOPLS. The other development branches are 
only available on the present GitLab repository.

# Method improvement

Once the method had been translated, verified and validated, it was tested 
on a real dataset.

(Coming soon)

The next step is to introduce variable selection to this ConsensusOPLS-DA 
method.

(Coming soon)
