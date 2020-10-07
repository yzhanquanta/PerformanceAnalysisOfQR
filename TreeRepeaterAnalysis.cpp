//
//  main.cpp
//  TreeRepeaterAnalysis
//
//  Created by YZhan on 9/1/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//

#include <iostream>
#include <math.h>
using namespace std;
#define testNum 10

double h(double x) {    // binary entropy function
    return -x*log2(x) - (1 - x)*log2(1 - x);
}

int main() {
    // explore possible tree structures for depth=3, 4, 5
    double etaD = 0.95; // detector effeciency
    // int numQR = 384;    // number of repeater nodes
    double C1Rec[testNum] = {0}; // store the calculated results
    double C2Rec[testNum] = {0};
    double R1Rec[testNum] = {0};
    double R2Rec[testNum] = {0};
    double L[testNum] = {0};    // distance between Alice and Bob
    double LStep = (2000.0 - 200.0)/(testNum - 1);
    for(int i = 0; i < testNum; ++i) {
        L[i] = 200 + LStep*i;
    }
    double Latt = 20;   // attenuation distance in fiber
    // double epR = 0.0003;    // error in re-encoding
    double epRArray[4] = {0};
    epRArray[0] = 0.0001;
    epRArray[1] = 0.0003;
    epRArray[2] = 0.0005;
    epRArray[3] = 0.001;
    for (int epRInd = 0; epRInd < 4; ++epRInd) {
        double epR = epRArray[epRInd];
        cout<< "epsilon_r = "<< epR<< endl;
        
        
        double tauPh = 1;  // time to generate a qubit, in units of ns
        double tauCZ = 10;  // spin-spin CZ gate time, in units of ns
        int level1TimeRatio = 100;  // time to generate level-1 photon = 100*tauPh
        
        int testRange = 15; // test range of branching parameters for each level
        int branPara[10] = {0};  // branching parameters
        double R[10] = {0}; // success prob. of indirect Z-measurement in level i (i<depth)
        
        for(int LInd = 0; LInd < testNum; ++LInd){
            // cout<< "L="<< L[LInd]<< " km, Borregaard's: ";
            double C1Opt = INFINITY; // best cost parameter
            double rs1Opt = 0;
            int numQR1Opt = 0;
            int branPara1Opt[10] = {0}; // best branching parameters
            int depth1Opt = 0; // best depth
            double C2Opt = INFINITY; // best cost parameter, 1 for PRX protocol, 2 for our protocol
            double rs2Opt = 0;
            int numQR2Opt = 0;
            int branPara2Opt[10] = {0}; // best branching parameters
            int depth2Opt = 0; // best depth
            for(int numQR = 1; numQR < L[LInd]; ++numQR) {
                // cout<< "Testing m = "<< numQR<< endl;
                // secret bit fraction
                double epTrans = 1 - pow((1 - epR), (numQR + 1));   // transmission error
                // double epTrans = (numQR + 1)*epR;
                double Q = 2*epTrans/3; // qubit error rate
                double f = 1 - h(Q) - Q - (1 - Q)*h((1 - 3*Q/2)/(1 - Q));
                // cout<< f<< endl;
                if(Q > 0.1261) {
                    // cout<< "The maximum number of repeater nodes is "<< numQR<< endl;
                    break;
                }
                double repeaterSep = L[LInd]/(numQR + 1); // distance between repeater nodes
                double eta = exp(-repeaterSep/Latt);    // transmission prob. between repeater nodes
                double mu = 1 - eta*etaD;   // single-photon loss rate
                for(int i = 0; i < testRange; ++i) {    // level-5, 0 means depth<=4
                    branPara[5] = i;    // b_4
                    R[5] = 0;
                    for(int j = 0; j < testRange; ++j) {    // level-4, 0 means depth<=3
                        branPara[4] = j;    // b_3
                        if((i > 0) && (j == 0)) { // b_4!=0, so b_3 cannot be 0
                            continue;;
                        }
                        R[4] = 1 - pow(1 - (1 - mu)*pow(1 - mu + mu*R[6], branPara[6]), branPara[5]);
                        for(int k = 1; k < testRange; ++k) {    // level-3
                            branPara[3] = k;    // b_2
                            R[3] = 1 - pow(1 - (1 - mu)*pow(1 - mu + mu*R[5], branPara[5]), branPara[4]);
                            for(int l = 1; l < testRange; ++l) {    // level-2
                                branPara[2] = l;    // b_1
                                R[2] = 1 - pow(1 - (1 - mu)*pow(1 - mu + mu*R[4], branPara[4]), branPara[3]);
                                for (int m = 1; m < testRange; ++m) {   // level-1
                                    branPara[1] = m;    // b_0
                                    R[1] = 1 - pow(1 - (1 - mu)*pow(1 - mu + mu*R[3], branPara[3]), branPara[2]);
                                    // calculate cost parameter C and secret key rate r_s
                                    // tree information
                                    int depth = 0;  // depth of the tree
                                    if(j == 0) {
                                        depth = 3;
                                    }
                                    else if(i == 0) {
                                        depth = 4;
                                    }
                                    else {
                                        depth = 5;
                                    }
                                    int totNumQ = 1;    // total # of qubits in the tree
                                    int nodeNum[10] = {0};   // number of nodes in each level
                                    nodeNum[0] = 1;
                                    for(int dd = 1; dd <= depth; ++dd) {
                                        nodeNum[dd] = nodeNum[dd - 1]*branPara[dd];
                                        totNumQ += nodeNum[dd];
                                    }
                                    // transmission probability
                                    double effLossRate = 1 - (pow(1 - mu + mu*R[1], branPara[1]) - pow(mu*R[1], branPara[1]))*pow(1 - mu + mu*R[2], branPara[2]);
                                    double transProb = pow((1 - effLossRate), (numQR + 1));
                                    // tree generation time
                                    double tGene1 = (totNumQ - 1 + (level1TimeRatio - 1)*branPara[1])*tauPh + (totNumQ - nodeNum[depth] - 1 + 2*branPara[1])*tauCZ;
                                    double tGene2 = ((nodeNum[depth] + nodeNum[depth - 1])*(depth + 1) - 2)*tauPh;
                                    // cout<< tGene1<< endl;
                                    // cost parameter and secret key rate
                                    double C1 = tGene1*numQR*(depth + 1)*Latt/(tauPh*L[LInd]*f*transProb);
                                    double rs1 = f*transProb/(tGene1*pow(10, -9));
                                    double C2 = tGene2*numQR*2*Latt/(tauPh*L[LInd]*f*transProb);
                                    double rs2 = f*transProb/(tGene2*pow(10, -9));
                                    if((C1 < C1Opt)) {  // for PRX protocol
                                        // if((C1 < C1Opt) && (totNumQ < 300)) {
                                        // if((rs1 > rs1Opt)) {
                                        C1Opt = C1;
                                        rs1Opt = rs1;
                                        depth1Opt = depth;
                                        numQR1Opt = numQR;
                                        for(int dd = 1; dd <= depth; ++dd) {
                                            branPara1Opt[dd] = branPara[dd];
                                        }
                                    }
                                    if((C2 < C2Opt)) {  // for our protocol
                                        // if((rs2 > rs2Opt)) {
                                        C2Opt = C2;
                                        rs2Opt = rs2;
                                        depth2Opt = depth;
                                        numQR2Opt = numQR;
                                        for(int dd = 1; dd <= depth; ++dd) {
                                            branPara2Opt[dd] = branPara[dd];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            C1Rec[LInd] = C1Opt;
            C2Rec[LInd] = C2Opt;
            R1Rec[LInd] = rs1Opt;
            R2Rec[LInd] = rs2Opt;
            // output
            //        // cout<< "For Borregaard's protocol, the optimized (minimum) cost parameter C = "<< C1Opt<< ", with a depth-"<< depth1Opt<< " tree, whose branching parameters are {";
            //        cout<< "C="<< C1Opt<< ", r_s="<< rs1Opt<< ", m="<< numQR1Opt<< ", with {";
            //        for(int dd = 1; dd < depth1Opt; ++dd) {
            //            cout<< branPara1Opt[dd]<< " ";
            //        }
            //        cout<< branPara1Opt[depth1Opt]<< "}. Ours: C="<< C2Opt<< ", r_s="<< rs2Opt<< ", m="<< numQR2Opt<< ", with {";
            //        // cout<< branPara1Opt[depth1Opt]<< "}, and the number of repeater nodes m = "<< numQR1Opt<< ", and the corresponding secret key rate r_s = "<< rs1Opt<< endl;
            //
            //        // cout<< "For our protocol, the optimized (minimum) cost parameter C = "<< C2Opt<< ", with a depth-"<< depth2Opt<< " tree, whose branching parameters are {";
            //        for(int dd = 1; dd < depth2Opt; ++dd) {
            //            cout<< branPara2Opt[dd]<< " ";
            //        }
            //        cout<< branPara2Opt[depth1Opt]<< "}"<< endl;
            //        // cout<< branPara2Opt[depth1Opt]<< "}, and the number of repeater nodes m = "<< numQR2Opt<< ", and the corresponding secret key rate r_s = "<< rs2Opt<< endl;
            cout<< L[LInd]<< " done;"<< endl;
        }
        
        cout<< "Borregaard's:"<< endl;
        cout<< "Cost parameters:"<< endl;
        for (int LInd = 0; LInd < testNum; ++LInd) {
            cout<< C1Rec[LInd]<< ",";
        }
        cout<< endl<< "Rate:"<< endl;
        for (int LInd = 0; LInd < testNum; ++LInd) {
            cout<< R1Rec[LInd]<< ",";
        }
        
        cout<< endl<< "Ours:"<< endl;
        cout<< "Cost parameters:"<< endl;
        for (int LInd = 0; LInd < testNum; ++LInd) {
            cout<< C2Rec[LInd]<< ",";
        }
        cout<< endl<< "Rate:"<< endl;
        for (int LInd = 0; LInd < testNum; ++LInd) {
            cout<< R2Rec[LInd]<< ",";
        }
    }
    return 0;
}
