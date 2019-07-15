#include <iostream>
#include <fstream>
#include <vector>
#include "matrix.hpp"

void infos(int const& c, int& s) {

}

int main() {
    std::ifstream ifile;
    ifile.open ("matrix.txt");
    std::vector<mtx> matrices;
    int h, w;
    while (!ifile.eof()) {
        ifile >> h >> w;
        mtx m(h, w);
        ifile >> m;
        matrices.push_back(m);
    }
    ifile.close();
    int i = 1;
    int c, s;
    while (1) {
        std::cout << "A futo program tartalmazott matrixai:" << std::endl << std::endl;
        for(auto it = matrices.begin(); it != matrices.end(); ++it)
            std::cout << i++ << "." << std::endl << *it << std::endl << std::endl << std::endl;

        std::cout << "A lehetseges muveletek:" << std::endl <<
                  '\t' << "1) osszeadas" << std::endl <<
                  '\t' << "2) matrixok szorzasa" << std::endl <<
                  '\t' << "3) skalarral szorzas" << std::endl <<
                  '\t' << "4) transzponalas" << std::endl <<
                  '\t' << "5) inverz matrix szamitasa" << std::endl <<
                  '\t' << "6) rang szamitasa" << std::endl << std::endl <<
                  '\t' << "7) determinans szamitasa" << std::endl <<
                  '\t' << "8) Gauss eliminacio elvegzese" << std::endl <<
                  "Egyeb:" << std::endl <<
                  '\t' << "9) uj matrix megadasa" << std::endl << std::endl << "Valasztas: ";

        std::cin >> c;
        i = 1;


        if (1 <= c && c <= 5) {
            std::cout << std::endl << "A muvelet az eredmenyt..." << std::endl <<
                      '\t' << "1) ne mentse el" << std::endl <<
                      '\t' << "2) uj matrixkent mentse el" << std::endl;
            if (3 <= c)
                std::cout << '\t' << "3) az eredeti matrixba mentse" << std::endl;
        }

        std::cout << std::endl << "Valasztas: ";

        std::cin >> s;

        std::cout << std::endl << "Valassz az alabbi matrixok kozul:" << std::endl;
        for(auto it = matrices.begin(); it != matrices.end(); ++it)
            std::cout << i++ << "." << std::endl << std::endl << *it << std::endl << std::endl;

        std::cout << std::endl << "Valasztas:" << std::endl;

        int a, b;

        if (1 <= c && c <= 8)
            std::cin >> a;

        if (1 == c || c == 2)
            std::cin >> b;
        --a;
        --b;
        mtx m;
        double D;
        std::cout << m.geth();
        switch(c) {
            case 1:
                m = matrices[a] + matrices[b];

                break;
            case 2:
                m = matrices[a] * matrices[b];
                break;
            case 3:
                std::cout << "Adj meg egy skalart: " << std::endl;
                double scalar;
                std::cin >> scalar;
                m = matrices[a] * scalar;
                break;
            case 4:
                m = matrices[a].transp();
                break;
            case 5:
                m = matrices[a].inverse();
                break;
            case 6:
                std::cout << matrices[a].r << std::endl;
                break;
            case 7:
                D = matrices[a].det();
                std::cout << D << std::endl;
                break;
            case 8:
                matrices[a].Gauss();
                break;
            default:
                break;
        }
        if (c == 9) {
            std::cout << std::endl << "Add meg a matrix dimenzioit (magassag, szelesseg sorrendben)" << std::endl;
            std::cin >> h >> w;
            mtx mx(h, w);
            std::cin >> mx;
            matrices.push_back(mx);
        }

        if (m.geth() != 0)
            std::cout << std::endl << "Az eredmeny:" << std::endl << m << std::endl;
        switch (s) {
            case 2:
                matrices.push_back(m);
                break;
            case 3:
                matrices[a] = m;
                break;
            default:
                break;
        }
        std::ofstream ofile;
        ofile.open ("matrix.txt");
        for(auto it = matrices.begin(); it != matrices.end(); ++it) {
            if (it != matrices.begin())
                ofile << std::endl << std::endl;
            ofile << (*it).geth() << std::endl << (*it).getw() << *it;
        }


        ofile.close();
    }
    return 0;
}
