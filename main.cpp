#define _CRT_RAND_S

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <random>
#include <cmath>
#include "gnuplot_iostream.h"
#include <numeric>
#include <cstdlib>
using namespace std;

const int m = 3;
const int n = 1;
const int M = 101;
const int fcm = 2;
const double eps = 1e-10;
const double iter = 10;

void updateD(vector<vector<double>>& D, vector<vector<double>>& spiralka,  vector<vector<double>>& V) {
    for (int j = 0; j < m; j++) {
        for (int k = 0; k < M; k++) {
            D[j][k] = sqrt(pow((spiralka[k][0] - V[j][0]),2) + pow((spiralka[k][1] - V[j][1]),2));
            if (D[j][k] < 1e-5)
                D[j][k] = 1e-5;
        }
    }
}

void updateU(vector<vector<double>>& D,vector<vector<double>>& U) {
    for (int k = 0; k < M; k++) {
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int l = 0; l < m; l++) {
                sum += pow( D[j][k] / D[l][k] , 2.0 / (fcm - 1) );
            }
            U[j][k] = 1.0 / sum;
        }
    }
}

void updateV(vector<vector<double>>& V,vector<vector<double>>& U,vector<vector<double>>& spiralka) {
    for (int x = 0; x < m; x++) {
        double sumU = 0;
        double sumX = 0;
        double sumY = 0;
        for (int s = 0; s < M; s++) {
            double u2 = pow(U[x][s], fcm);
            sumU += u2;
            sumX += u2 * spiralka[s][0]; // x
            sumY += u2 * spiralka[s][1]; // y
        }

        V[x][0] = sumX / sumU;
        V[x][1] = sumY / sumU;
    }
}

void coloring(vector<vector<double>> U, vector<vector<double>>& spiralka) {
    for (int k = 0; k < M; k++) {
        double h = 0.0;
        for (int j = 0; j < m; j++) {
            h += U[j][k] * j;
        }
        spiralka[k][2] = h;
    }
}

int main() {
    unsigned int number;
    vector<vector<double>> spiralka;
    string line = "", first = "", second = ""; // reading csv rariables
    int location = 0, readcount = 0;
    Gnuplot gp("D:/gnuplot/bin/gnuplot.exe");
    ifstream iFile("../spiralka.csv");

    while (getline(iFile, line) && readcount < M)
    {
        double group = readcount%4; // temporary colour
        double group_odl = readcount%4; // temporary colour
        location = line.find(','); // spliting line with,
        first = line.substr(0, location);
        line = line.substr(location + 1, line.length());
        spiralka.push_back({stof(first), stof(line),group});
        readcount++;
    }

    vector<vector<double>> U;
    vector<vector<double>> D;
    vector<vector<double>> V;

    for (int i = 0; i < m; i++) {
        vector<double> rnd;
        for (int j = 0; j < M; j++) {
            rand_s(&number);
            rnd.push_back(static_cast<double>(number) / ((static_cast<double>(UINT_MAX) + 1) * 1.1) + 0.1);
        }
        D.push_back(rnd);
    }

    for (int i = 0; i < m; i++) {
        double resD = 0;
        vector<double> resDik;
        for (int j = 0; j < M; j++) {
        resD += pow(D[i][j], 1/(1 - fcm));
        }
        for (int j = 0; j < M; j++) {
            resDik.push_back(pow(D[i][j], 1/(1 - fcm)) / resD);
        }
        U.push_back(resDik);
    }

    for (int i = 0; i < m; i++) {
        double sumU = 0;
        double sumX = 0;
        double sumY = 0;
        for (int s = 0; s < M; s++) {
            double u2 = pow(U[i][s], fcm);
            sumU += u2;
            sumX += u2 * spiralka[s][0]; // x
            sumY += u2 * spiralka[s][1]; // y
        }
        V.push_back({sumX/sumU, sumY/sumU, static_cast<double>(i)});
    }

    for (int i = 0; i < iter; i++) {
        updateD(D,spiralka,V);
        updateU(D,U);
        updateV(V,U,spiralka);
    }

    coloring(U, spiralka);

    gp << "set term wxt 1\n";
    gp << "set title 'Wykres spiralka z centroidami'\n";
    gp << "set cbrange [0:2]\n";
    gp << "set palette rgb 3,13,10\n";
    gp << "plot '-' using 1:2:3 with points pt 7 lc palette title 'punkty', "
             "'-' using 1:2:3 with points pt 6 ps 1.5 lc palette title 'środki'\n";

    gp.send1d(spiralka);
    gp.send1d(V);
    cin.get();
}