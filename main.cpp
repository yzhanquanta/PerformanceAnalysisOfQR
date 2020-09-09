//
//  main.cpp
//  EffErrorProb
//
//  Created by YZhan on 9/8/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//

#include <iostream>
#include <math.h>
using namespace std;

class tree {
private:
    int depth = 0;  // depth of the tree, assumed to be [2, 10]
    int branPara[10] = {0}; // branching parameters
    
public:
    tree() {}   // consturctor with no variable
    tree(int d, int *b) {   // constructor with depth and branching parameters
        if(d > 10) {
            cout<< "Depth out of reach!"<< endl;
        }
        depth = d;
        for(int dd = 0; dd < depth; ++dd) {
            if(b[dd] > 24) {
                cout<< "Branching parameter out of reach!"<< endl;
            }
            branPara[dd] = b[dd];
        }
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
    
    ~tree() {}
    
};

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

int main() {
    // tree information
    int d1 = 3;  // assume depth >= 2
    int b1[10] = {4, 14, 4, 0, 0, 0, 0, 0, 0, 0};
    tree tree1(d1, b1);
    
    // effective loss prob and indirect meas prob
    double Latt = 20;   // attenuation distance (km)
    double L = 1000;    // distance between Alice and Bob (km)
    int numQR = 402;  // number of repeater nodes
    double etaD = 0.95; // detector efficiency
    double mu = 1 - exp(-L/(numQR + 1)/Latt)*etaD;  // single-photon loss rate
    // double mu = 0.2;
    
    double R[10] = {0}; // indirect Z success prob
    double s[10] = {0}; // single indirect Z success prob
    for(int dd = tree1.getDepth() - 1; dd >= 0; --dd) {
        s[dd] = (1 - mu)*pow(1 - mu + mu*R[dd + 2], tree1.getBranPara(dd + 1));
        R[dd] = 1 - pow(1 - s[dd], tree1.getBranPara(dd));
        // cout<< dd<< " "<< s[dd]<< endl;
    }
    
    // recursive calculations of e_{I_k}, e_{I_k|m_k}, p_k(m_k), e_{I_k|1}, and e_{n_k}
    double P_indir_Z[10] = {0}; // conditional prob of an indirect Z meas given a success of Z meas
    for(int dd = tree1.getDepth() - 1; dd >= 0; --dd) {
        P_indir_Z[dd] = R[dd]/(1 - mu + mu*R[dd]);
        // cout<< dd<< " "<< P_indir_Z[dd]<< endl;
    }
    double epsilon = 1.0*pow(10, -4);    // single-photon error prob
    cout<< "The single-photon error prob. epsilon = "<< epsilon<< endl;
    double e_Ik[10] = {0};  // given by Eq.13, the error prob of an indirect meas on level-k
    double e_Ik_mk[10][25] = {0};   // given by Eq.15, the error prob of an indirect meas on level-k given m_k branches measured, assuming max=25=testRange
    double p_k_mk[10][25] = {0};    // given by Eq.14
    double e_Ik_1[10] = {0};    // given by Eq.16, the error prob of an indirect meas on level-k given by a single branch
    double e_nk[10][25] = {0};  // given by Eq.17, the error prob of an indirect meas on level-k given by a single branch where n_k level-(k+2) qubits are directly measured
    
    e_Ik_1[tree1.getDepth() - 1] = epsilon;    // boundary conditions are d-1 and d-2, calculate e_{I_k|1} for d-1 and d-2
    for(int i = 0; i <=  tree1.getBranPara(tree1.getDepth() - 1) + 1; ++i) {
        if((i % 2) == 0)
            continue;;
        e_Ik_1[tree1.getDepth() - 2] += permutation(tree1.getBranPara(tree1.getDepth() - 1) + 1, i)*pow(epsilon, i)*pow(1 - epsilon, tree1.getBranPara(tree1.getDepth() - 1) + 1 - i);
    }
    
    for(int dd = tree1.getDepth() - 2; dd < tree1.getDepth(); ++dd) {
        // calculate e_I_k_m_k for d-1 and d-2
        for(int mk = 1; mk <= tree1.getBranPara(dd); ++mk) {
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
        for(int mk = 1; mk <= tree1.getBranPara(dd); ++mk) {
            p_k_mk[dd][mk] = permutation(tree1.getBranPara(dd), mk)*pow(s[dd], mk)*pow(1 - s[dd], tree1.getBranPara(dd) - mk);
        }
        // calculate e_Ik for d-1 and d-2
        for(int mk = 1; mk <= tree1.getBranPara(dd); ++mk) {
            e_Ik[dd] += p_k_mk[dd][mk]*e_Ik_mk[dd][mk];
        }
        e_Ik[dd] = e_Ik[dd]/R[dd];
    }
    
    for(int dd = tree1.getDepth() - 3; dd >= 0; --dd) {    // recursively calculate all upper levels
        // calculate e_nk for level-dd
        for(int nk = 0; nk <= tree1.getBranPara(dd + 1); ++nk) {
            for(int i = 0; i <= nk + 1; ++i) {
                double temp = 0;
                for(int j = 0; j <= tree1.getBranPara(dd + 1) - nk; ++j) {
                    if((i + j) % 2 == 0)    // parity
                        continue;
                    temp += permutation(tree1.getBranPara(dd + 1) - nk, j)*pow(e_Ik[dd + 2], j)*pow(1 - e_Ik[dd + 2], tree1.getBranPara(dd + 1) - nk - j);
                }
                e_nk[dd][nk] += permutation(nk + 1, i)*pow(epsilon, i)*pow(1 - epsilon, nk + 1 - i)*temp;
            }
        }
        // calculate e_Ik_1 for level-dd
        for(int nk = 0; nk <= tree1.getBranPara(dd + 1); ++nk) {
            e_Ik_1[dd] += permutation(tree1.getBranPara(dd + 1), nk)*pow(P_indir_Z[dd + 2], tree1.getBranPara(dd + 1) - nk)*pow(1 - P_indir_Z[dd + 2], nk)*e_nk[dd][nk];
        }
        // calculate e_Ik_mk for level-dd
        for(int mk = 1; mk <= tree1.getBranPara(dd); ++mk) {
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
        for(int mk = 1; mk <= tree1.getBranPara(dd); ++mk) {
            p_k_mk[dd][mk] = permutation(tree1.getBranPara(dd), mk)*pow(s[dd], mk)*pow(1 - s[dd], tree1.getBranPara(dd) - mk);
        }
        // calculate e_Ik for level-dd
        for(int mk = 1; mk <= tree1.getBranPara(dd); ++mk) {
            e_Ik[dd] += p_k_mk[dd][mk]*e_Ik_mk[dd][mk];
        }
        e_Ik[dd] = e_Ik[dd]/R[dd];
    }
    
    // effective error prob of Z meas in tree-encoded RGS
    double e_n[25] = {0};   // given by Eq.12
    double e_Z = 0; // effective error prob of Z meas of a logical qubit
    for(int n = 0; n <= tree1.getBranPara(0); ++n) {
        // calculate e_n for Z meas of a logical qubit in RGS
        for(int i = 0; i <= n; ++i) {
            double temp = 0;
            for(int j = 0; j <= tree1.getBranPara(0) - n; ++j) {
                if((i + j) % 2 == 0)    // parity
                    continue;
                temp += permutation(tree1.getBranPara(0) - n, j)*pow(e_Ik[1], j)*pow(1 - e_Ik[1], tree1.getBranPara(0) - n - j);
            }
            e_n[n] += permutation(n, i)*pow(epsilon, i)*pow(1 - epsilon, n - i)*temp;
        }
    }
    for(int n = 0; n <= tree1.getBranPara(0); ++n) {
        e_Z += permutation(tree1.getBranPara(0), n)*pow(P_indir_Z[1], tree1.getBranPara(0) - n)*pow(1 - P_indir_Z[1], n)*e_n[n];
    }
    cout<< "RGS: The effective error prob. of Z-measurement on the logical qubit is e_Z = "<< e_Z<< endl;
    
    // effective error prob of X meas in tree-encoded RGS
    double e_X = e_Ik[0];
    cout<< "RGS: The effective error prob. of X-measurement on the logical qubit is e_X = "<< e_X<< endl;
    
    // effective error prob of tree-repeater and MBQC
    double P_succ = (pow(1 - mu + mu*R[1], tree1.getBranPara(0)) - pow(mu*R[1], tree1.getBranPara(0)))*pow(1 - mu + mu*R[2], tree1.getBranPara(1));    // success prob
    double P_n[25] = {0};   // given by Eq.24
    double P_error_n[25] = {0}; // given by Eq.25
    double e_decoding = 0;  // given by Eq.26, effective error prob
    double e_m[25] = {0};   // given by Eq.12
    double e_n2[25] = {0};  // given by Eq.25
    
    for(int m = 0; m <= tree1.getBranPara(1); ++m) {
        // calculate e_m
        for(int i = 0; i <= m; ++i) {
            double temp = 0;
            for(int j = 0; j <= tree1.getBranPara(1) - m; ++j) {
                if((i + j) % 2 == 0)    // parity
                    continue;
                temp += permutation(tree1.getBranPara(1) - m, j)*pow(e_Ik[2], j)*pow(1 - e_Ik[2], tree1.getBranPara(1) - m - j);
            }
            e_m[m] += permutation(m, i)*pow(epsilon, i)*pow(1 - epsilon, m - i)*temp;
        }
    }
    
    for(int n = 1; n <= tree1.getBranPara(0); ++n) {
        // calculate e_n2
        for(int i = 0; i <= n - 1; ++i) {
            double temp = 0;
            for(int j = 0; j <= tree1.getBranPara(0) - n; ++j) {
                if((i + j) % 2 == 0)    // parity
                    continue;
                temp += permutation(tree1.getBranPara(0) - n, j)*pow(e_Ik[1], j)*pow(1 - e_Ik[1], tree1.getBranPara(0) - n - j);
            }
            e_n2[n] += permutation(n - 1, i)*pow(epsilon, i)*pow(1 - epsilon, n - 1 - i)*temp;
        }
    }
    
    for(int n = 1; n <= tree1.getBranPara(0); ++n) {
        // calculate P_n
        P_n[n] = permutation(tree1.getBranPara(0), n)*pow(1 - mu, n)*pow(mu*R[1], tree1.getBranPara(0) - n)*pow(1 - mu + mu*R[2], tree1.getBranPara(1));
        // calculate P_error_n
        double e_Z_b1 = 0;  // given by Eq.25
        for(int m = 0; m <= tree1.getBranPara(1); ++m) {
            e_Z_b1 += permutation(tree1.getBranPara(1), m)*pow(P_indir_Z[2], tree1.getBranPara(1) - m)*pow(1 - P_indir_Z[2], m)*e_m[m];
        }
        P_error_n[n] = 1 - (1 - epsilon)*(1 - e_n2[n])*(1 - e_Z_b1);
    }
    double P_wrong = 0; // given by Eq.23
    for(int n = 1; n <= tree1.getBranPara(0); ++n) {
        P_wrong += P_n[n]*P_error_n[n];
    }
    e_decoding = P_wrong/P_succ;
    cout<< "Tree-repeater & MBQC: The effective error prob. is e_decoding = "<< e_decoding<< endl;
    
    // test
//    double test_e_I1_1 = 0;
//    for(int i = 0; i <= 3; ++i) {
//        test_e_I1_1 += permutation(3, i)*pow(epsilon, i)*pow(1 - epsilon, 3 - i);
//    }
//    cout<< test_e_I1_1<< endl;
    
    return 0;
}
