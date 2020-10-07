//
//  main.cpp
//  RGSPerf
//
//  Created by YZhan on 9/30/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//

#include <iostream>
#include <math.h>
using namespace std;
#define testNum 15

int factorial(int n) {
    int temp = 1;
    for(int i = 1; i <= n; ++i) {
        temp = temp*i;
    }
    return temp;
}

int permutation(int n, int m) {
    if(n < m) {
        cout<< "Incorrect input for number of permutation!"<< endl;
    }
    return factorial(n)/factorial(m)/factorial(n - m);
}

double binaryEntropy(double F) {
    return -F*log2(F) - (1 - F)*log2(1 - F);
}

class RGS {
private:
    int numBell = 0;    // number of branches of RGS N=2*numBell
    // encoding tree info
    int numMatterPerNode = 0;  // number of matter qubits per repeater
    int depth = 0;  // depth of the tree
    int branPara[10] = {0}; // branching parameters
    int maxBranPara = 20;   // maximum possible branching parameter
    int numPhotonLevel[10] = {0};   // number of photons at each level
    double mu = 0;  // single-photon loss rate
    double s[10] = {0}; // success prob of a single indirect Z-meas
    double R[10] = {0}; // success prob of indirect Z-meas
    double geneTime = 0;    // generation time of the tree (in units of Tcz)
    // logic meas
    double probLogicX = 0;  // success prob of a logic X meas
    double probLogicZ = 0;  // success prob of a logic Z meas
    double probBell = 0;    // success prob of a Bell meas
    // effective error prob of logic X and Z meas
    double logicErrorX = 0;
    double logicErrorZ = 0;
    // fidelity of entangled pair
    double fAB = 0;
    
public:
    RGS() {}   // consturctor with no variable
    
    RGS(int m, int d, int *b, double lossRate) {   // constructor with depth and branching parameters
        numBell = m;
        if(d > 10) {
            cout<< "Depth out of reach!"<< endl;
        }
        depth = d;
        mu = lossRate;
        
        // store branching parameters
        numPhotonLevel[0] = 1;
        for(int dd = 0; dd < depth; ++dd) {
            if(b[dd] > maxBranPara) {
                cout<< "Branching parameter out of reach!"<< endl;
            }
            branPara[dd] = b[dd];
            numPhotonLevel[dd + 1] = numPhotonLevel[dd]*branPara[dd];
        }
        
        // calculate s's and R's
        for(int dd = depth - 1; dd >= 0; --dd) {
            s[dd] = (1 - mu)*pow(1 - mu + mu*R[dd + 2], branPara[dd + 1]);
            R[dd] = 1 - pow(1 - s[dd], branPara[dd]);
        }
        
        // calculate the generation time of Buterako's
//        numMatterPerNode = depth + 1;
//        for (int dd = 1; dd < depth; ++dd) {
//            geneTime += numPhotonLevel[dd];
//        }
//        geneTime = 2*numBell*(2 + geneTime);
        
        // calculate the generation time of ours
        numMatterPerNode = 1;
        geneTime = 2*numBell*(numPhotonLevel[depth] + numPhotonLevel[depth - 1])*(depth + 3) - 2;
        
        // calculate success prob of logic X and Z meas
        probLogicX = R[0];
        probLogicZ = pow(1 - mu + mu*R[1], branPara[0]);
        probBell = pow(1 - mu, 2)/2;
    }
    
    void errorAnalysis(double epsilonX, double epsilonZ) {
        if ((epsilonX == 0) && (epsilonZ == 0)) {   // the X and Z meas prob have been given in constructor
            return;
        }
        
        // calculate the success prob in presence of single-photon error
        double P_indir_Z[10] = {0}; // conditional prob of an indirect Z meas given a success of Z meas
        for(int dd = depth - 1; dd >= 0; --dd) {
            P_indir_Z[dd] = R[dd]/(1 - mu + mu*R[dd]);
            // cout<< dd<< " "<< P_indir_Z[dd]<< endl;
        }
        // cout<< "The single-photon X-meas error prob epsilonX = "<< epsilonX<< ", Z-meas error prob epsilonZ = "<< epsilonZ<< endl;
        double e_Ik[10] = {0};  // given by Eq.13, the error prob of an indirect meas on level-k
        double e_Ik_mk[10][35] = {0};   // given by Eq.15, the error prob of an indirect meas on level-k given m_k branches measured, assuming max=35=testRange
        double p_k_mk[10][35] = {0};    // given by Eq.14
        double e_Ik_1[10] = {0};    // given by Eq.16, the error prob of an indirect meas on level-k given by a single branch
        double e_nk[10][35] = {0};  // given by Eq.17, the error prob of an indirect meas on level-k given by a single branch where n_k level-(k+2) qubits are directly measured
        
        e_Ik_1[depth - 1] = epsilonX;    // boundary conditions are d-1 and d-2, calculate e_{I_k|1} for d-1 and d-2
        for(int i = 0; i <=  branPara[depth - 1]; ++i) {    // i is the number of Z errors in level-d
            if((i % 2) == 0) {    // if even, X meas is wrong
                e_Ik_1[depth - 2] += epsilonX*permutation(branPara[depth - 1], i)*pow(epsilonZ, i)*pow(1 - epsilonZ, branPara[depth - 1] - i);
            }
            else {  // if odd, only Z meas errors
                e_Ik_1[depth - 2] += (1 - epsilonX)*permutation(branPara[depth - 1], i)*pow(epsilonZ, i)*pow(1 - epsilonZ, branPara[depth - 1] - i);
            }
        }
        
        for(int dd = depth - 2; dd < depth; ++dd) {
            // calculate e_I_k_m_k for d-1 and d-2
            for(int mk = 1; mk <= branPara[dd]; ++mk) {
                if((mk % 2) == 1) {   // mk is odd
                    for(int i = (mk + 1)/2; i <= mk; ++i) { // mojority vote
                        e_Ik_mk[dd][mk] += permutation(mk, i)*pow(e_Ik_1[dd], i)*pow(1 - e_Ik_1[dd], mk - i);
                    }
                }
                else {  // mk is even
                    for(int i = mk/2; i <= mk - 1; ++i) {
                        e_Ik_mk[dd][mk] += permutation(mk - 1, i)*pow(e_Ik_1[dd], i)*pow(1 - e_Ik_1[dd], mk - 1 - i);
                    }
                }
            }
            // calculate p_k(m_k) for d-1 and d-2
            for(int mk = 1; mk <= branPara[dd]; ++mk) {
                p_k_mk[dd][mk] = permutation(branPara[dd], mk)*pow(s[dd], mk)*pow(1 - s[dd], branPara[dd] - mk);
            }
            // calculate e_Ik for d-1 and d-2
            for(int mk = 1; mk <= branPara[dd]; ++mk) {
                e_Ik[dd] += p_k_mk[dd][mk]*e_Ik_mk[dd][mk];
            }
            e_Ik[dd] = e_Ik[dd]/R[dd];
        }
        
        for(int dd = depth - 3; dd >= 0; --dd) {    // recursively calculate all upper levels
            // calculate e_nk for level-dd
            for(int nk = 0; nk <= branPara[dd + 1]; ++nk) {
                for(int i = 0; i <= nk; ++i) {
                    double temp = 0;
                    for(int j = 0; j <= branPara[dd + 1] - nk; ++j) {
                        if((i + j) % 2 == 0) {    // parity
                            temp += epsilonX*permutation(nk, i)*pow(epsilonZ, i)*pow(1 - epsilonZ, nk - i)*permutation(branPara[dd + 1] - nk, j)*pow(e_Ik[dd + 2], j)*pow(1 - e_Ik[dd + 2], branPara[dd + 1] - nk - j);
                        }
                        else {
                            temp += (1 - epsilonX)*permutation(nk, i)*pow(epsilonZ, i)*pow(1 - epsilonZ, nk - i)*permutation(branPara[dd + 1] - nk, j)*pow(e_Ik[dd + 2], j)*pow(1 - e_Ik[dd + 2], branPara[dd + 1] - nk - j);
                        }
                    }
                    e_nk[dd][nk] += temp;
                }
            }
            // calculate e_Ik_1 for level-dd
            for(int nk = 0; nk <= branPara[dd + 1]; ++nk) {
                e_Ik_1[dd] += permutation(branPara[dd + 1], nk)*pow(P_indir_Z[dd + 2], branPara[dd + 1] - nk)*pow(1 - P_indir_Z[dd + 2], nk)*e_nk[dd][nk];
            }
            // calculate e_Ik_mk for level-dd
            for(int mk = 1; mk <= branPara[dd]; ++mk) {
                if(mk % 2 == 1) {   // mk is odd
                    for(int i = (mk + 1)/2; i <= mk; ++i) { // mojority vote
                        e_Ik_mk[dd][mk] += permutation(mk, i)*pow(e_Ik_1[dd], i)*pow(1 - e_Ik_1[dd], mk - i);
                    }
                }
                else {  // mk is even
                    for(int i = mk/2; i <= mk - 1; ++i) {
                        e_Ik_mk[dd][mk] += permutation(mk - 1, i)*pow(e_Ik_1[dd], i)*pow(1 - e_Ik_1[dd], mk - 1 - i);
                    }
                }
            }
            // calculate p_k(m_k) for level-dd
            for(int mk = 1; mk <= branPara[dd]; ++mk) {
                p_k_mk[dd][mk] = permutation(branPara[dd], mk)*pow(s[dd], mk)*pow(1 - s[dd], branPara[dd] - mk);
            }
            // calculate e_Ik for level-dd
            for(int mk = 1; mk <= branPara[dd]; ++mk) {
                e_Ik[dd] += p_k_mk[dd][mk]*e_Ik_mk[dd][mk];
            }
            e_Ik[dd] = e_Ik[dd]/R[dd];
        }
        
        // effective error prob of Z meas in tree-encoded RGS
        double e_n[35] = {0};   // given by Eq.12
        double e_Z = 0; // effective error prob of Z meas of a logical qubit
        for(int n = 0; n <= branPara[0]; ++n) {
            // calculate e_n for Z meas of a logical qubit in RGS
            for(int i = 0; i <= n; ++i) {
                double temp = 0;
                for(int j = 0; j <= branPara[0] - n; ++j) {
                    if((i + j) % 2 == 0)    // parity
                        continue;
                    temp += permutation(branPara[0] - n, j)*pow(e_Ik[1], j)*pow(1 - e_Ik[1], branPara[0] - n - j);
                }
                e_n[n] += permutation(n, i)*pow(epsilonZ, i)*pow(1 - epsilonZ, n - i)*temp;
            }
        }
        for(int n = 0; n <= branPara[0]; ++n) {
            e_Z += permutation(branPara[0], n)*pow(P_indir_Z[1], branPara[0] - n)*pow(1 - P_indir_Z[1], n)*e_n[n];
        }
        // cout<< "RGS: The effective error prob. of Z-measurement on the logical qubit is e_Z = "<< e_Z<< endl;
        
        // effective error prob of X meas in tree-encoded RGS
        double e_X = e_Ik[0];
        // cout<< "RGS: The effective error prob. of X-measurement on the logical qubit is e_X = "<< e_X<< endl;
        
        logicErrorX = e_X;
        logicErrorZ = e_Z;
        
        // test
        //        for(int dd = 0; dd < depth; ++dd) {
        //            cout<< e_Ik[dd]<< endl;
        //        }
    }
    
    double fidelity(double epsilonX, double epsilonZ, int numQR) {
        errorAnalysis(epsilonX, epsilonZ);
        double Ex = 0.25*(1 - pow(1 - 2*epsilonX, 2*(numQR + 1))*pow(1 - 2*logicErrorX, 2*numQR));
        double Ez = Ex;
        double Ey = 0.25*(1 + pow(1 - 2*epsilonX, 2*(numQR + 1))*pow(1 - 2*logicErrorX, 2*numQR) - 2*pow(1 - 2*epsilonX, 2*(numQR + 1))*pow(1 - 2*logicErrorX, numQR)*pow(1 - 2*logicErrorZ, (2*numBell - 2)*numQR));
        
        return 1 - (Ex + Ey + Ez);
    }
    
    double getErrorX() {
        return logicErrorX;
    }
    
    double getErrorZ() {
        return logicErrorZ;
    }
    
    double probSucc(int numQR) {
        return pow((1 - pow(1 - probBell, numBell))*pow(probLogicX, 2)*pow(probLogicZ, 2*numBell - 2), numQR + 1);
    }
    
    double rateRGS(int numQR) { // the communication rate of RGS (in units of 1/Tcz)
        double pRGS = (1 - pow(1 - probBell, numBell))*pow(probLogicX, 2)*pow(probLogicZ, 2*numBell - 2);  // success prob of connecting two RGS
        return pow(pRGS, numQR + 1)/geneTime;
    }
    
    double rateRGSPerMatter(int numQR) {
        return rateRGS(numQR)/(numQR*numMatterPerNode);
    }
    
    double rateSKPerMatter(double epsilonX, double epsilonZ, int numQR) {
        return rateRGSPerMatter(numQR)*(1 - 2*binaryEntropy(fidelity(epsilonX, epsilonZ, numQR)));
    }
    
    int getNumBell() {
        return numBell;
    }
    
    int getDepth() {    // get the depth
        return depth;
    }
    
    int getBranPara(int level) {    // get a branching parameter
        return branPara[level];
    }
    
    void print() {
        cout<< "The depth of this tree is "<< depth<< ", with branching parameters {";
        for(int dd = 0; dd < depth - 1; ++dd) {
            cout<< branPara[dd]<< ", ";
        }
        cout<< branPara[depth - 1]<< "}."<< endl;
    }
    
    ~RGS() {}
};

int main() {
    // basic information
    double L[testNum] = {0};
    double rateOpt[testNum] = {0};
    double LStep = (200.0 - 10.0)/(testNum - 1);
    for(int i = 0; i < testNum; ++i) {
        L[i] = 10 + LStep*i;
    }
    // double L = 10;    // distance between Alice and Bob (in units of Latt=20 km)
    double Latt = 20;
    double etaCD = 1;   // photon collection and detector efficiency \eta_c\eta_d
    double epsilon = 0.00001;    // single-photon error prob
    
    // single-case test
    //    int d = 2;  // depth of the tree
    //    int b[10] = {7, 4, 0, 0, 0, 0, 0, 0, 0, 0};    // branching parameters of the tree
    //    int numQR = 49;    // number of repeater stations in between
    //    int numBell = 11;   // number of Bell pairs in a meas node (N=2*numBell)
    //    double mu = 1 - etaCD*exp(-L/(numQR + 1)/2);  // single-photon loss rate
    //    RGS testRGS(numBell, d, b, mu);
    //    testRGS.errorAnalysis(0.001, 0.001);
    //    cout<< "RGS rate: "<< testRGS.rateRGSPerMatter(numQR)<< endl;
    //    cout<< "Error probability of logic X and Z measurements: X - "<< testRGS.getErrorX()<< ", Z - "<< testRGS.getErrorZ()<< endl;
    
    
    for(int LInd = 0; LInd < testNum; ++LInd) {
        // optimum searching
        double optimalRate = -INFINITY;
        int optimalNumBell = 0;
        int optimalNumQR = 0;
        int optimalBranPara[10] = {0};
        int optimalDepth = 0;
        
        int testRange = 11; // test range for each level
        int testRangeBell = 20; // test range of number of Bell pairs at each meas node
        int b[10] = {0}; // branching parameters
        
        // explore depth = 2-3 to final the optimum
        for(int l = 1; l < testRange; ++l) {    // level-2
            b[1] = l;
            for (int m = 1; m < testRange; ++m) {   // level-1
                b[0] = m;
                
                // tree depth
                int d = 2;  // depth of the tree
                
                // find optimal numQR and numBell
                for (int numQR = 1; numQR < L[LInd]*Latt; ++numQR) {
                    double mu = 1 - etaCD*exp(-L[LInd]/(numQR + 1)/2);  // single-photon loss rate
                    for (int numBell = 2; numBell < testRangeBell; ++numBell) {
                        RGS testRGS(numBell, d, b, mu);
                        // double testRate = testRGS.rateRGSPerMatter(numQR);
                        double testRate = testRGS.rateSKPerMatter(epsilon, epsilon, numQR);
                        if (testRate > optimalRate) {
                            optimalRate = testRate;
                            optimalDepth = d;
                            optimalNumQR = numQR;
                            optimalNumBell = numBell;
                            for (int dd = 0; dd < optimalDepth; ++dd) {
                                optimalBranPara[dd] = b[dd];
                            }
                        }
                    }
                }
                
            }
        }
        
        rateOpt[LInd] = optimalRate;
        
        // output the result
        cout<< "For L = "<< L[LInd]*Latt<< " km,"<< endl;
        cout<< "The optimal RGS rate: "<< optimalRate<< endl;
        // cout<< "The optimal secret key rate: "<< optimalRate<< endl;
        cout<< "The optimal number of Bell pairs: "<< optimalNumBell<< endl;
        cout<< "The optimal number of repeater nodes: "<< optimalNumQR<< ", i.e. distance between repeaters: "<< L[LInd]/(optimalNumQR + 1)<< endl;
        cout<< "The optimal tree structure: {";
        for (int dd = 0; dd < optimalDepth - 1; ++dd) {
            cout<< optimalBranPara[dd]<< " ";
        }
        cout<< optimalBranPara[optimalDepth - 1]<< "}"<< endl;
        
    }
    
    for(int LInd = 0; LInd < testNum; ++LInd) {
        cout<< rateOpt[LInd]<< ",";
    }
    cout<< endl;
    
    return 0;
}
