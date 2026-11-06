#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <unistd.h>
#include <fstream>

using namespace std;

namespace Physics {
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    constexpr long double h = 6.62607015e-34L;
    constexpr long double hbar = h / (2.0L * pi);
}

class StationaryState{
    private: 
        const int n;
        const double E;
        vector<double> psi_spatial;

        double CalculateEnergy(int n, double L, double m) {
            double E = (Physics::hbar * n * Physics::pi) * (Physics::hbar * n * Physics::pi) / (2 * m * L * L);
            return E;
        }

        vector<double> CalculatePsiSpatial(int n, double L, int N) {
            vector<double> psi;
            psi.reserve(N);

            double dx = L / static_cast<double>(N);

            for(int i = 0; i < N; ++i) {
                double x = i * dx;
                psi.push_back(sqrt(2.0/L) * sin(n * Physics::pi / L * x));
            }
            return psi;
        }
    public:

        StationaryState(int n, double L, int N, double m) :
            n(n),
            E(CalculateEnergy(n, L, m)),
            psi_spatial(CalculatePsiSpatial(n, L, N)) {}
        

        int getN() const {
            return n;
        }

        double getEnergy() const {
            return E;
        }

        const vector<double>& getPsi() const {
            return psi_spatial;
        }
};

struct StateParams {
    double phi;
    double c;
    inline static double acc = 0.0;

    StateParams(double phi, double c) : phi(phi), c(c) {
        acc += c*c;
    }

    static void normalizeVector(vector<StateParams>& states) {
        if (acc == 0.0) return;

        for (auto& state : states) {
            state.c /= sqrt(acc);
        }
    }

    static void resetAcc() {
        acc = 0.0;
    }
};

class QuantumSystem{
    private:
        vector<StationaryState> states;
        vector<StateParams> params;
        
        const double L;
        const int intervals;
        const double mass;
    public:
        QuantumSystem(double L, int intervals, double mass):
            L(L), intervals(intervals), mass(mass) {}

        void addState(int n, double phi, double c) {
            states.emplace_back(n, L, intervals, mass);
            params.emplace_back(phi, c);
        }
        void prepareSystem() {
            StateParams::normalizeVector(params);
        }
        void simulate(double dt, double t_max) {
            cout << "Starting the simulation...\n";
            ofstream file("simulation.csv");
            if (!file.is_open()) return;

            vector<complex<double>> psiTotal(intervals);
            int step = 0;
            

            for (double t = 0.0; t <= t_max; t += dt) {
                fill(psiTotal.begin(), psiTotal.end(), complex<double>(0.0, 0.0));
                
                for(size_t k = 0; k < states.size(); ++k) {
                    const auto& psi_spatial = states[k].getPsi();
                    double E = states[k].getEnergy();
                    double c = params[k].c;
                    double phi = params[k].phi;

                    double omega = E / Physics::hbar * t;
                    complex<double> timePhase = polar(c, -(omega - phi));

                    for(int x = 0; x < intervals; x++) {
                        psiTotal[x] += timePhase * psi_spatial[x];
                    }

                }
                if(step % 10 == 0) {
                    file << t;
                    for (int i = 0; i < intervals; ++i) {
                        double probability_density = norm(psiTotal[i]); 
                        file << "," << probability_density;
                    }
                    file << "\n";
                    // cout << t << endl;
                }
                step++;

            }
            

            file.close();
            cout << "Simulation finished successfully\n";
        }

        void printParams() {
            cout << "=-=-=-= Quantum System Parameters =-=-=-=\n";
            cout << "Well lenght: " << L <<endl;
            cout << "Particle mass: " << mass << endl;

            for (auto& param : params) {
                cout << "|c|: "<<param.c << " phi: " << param.phi<<endl;
            }
        }
};

int enterN() {
    int N;
    cout << "Set the number of states: ";
    cin >> N;
    return N;
}


void Usage(char* argv[]) {
    cout << "Usage: " << argv[0] << " [L] [N] [dt] [tMax]\n";
}

int main(int argc, char* argv[]) {
    double L = 10.0;
    int intervals = 500;
    double dt = 0.01;
    double mass = 9.110e-31;
    double tMax = 10.0;

    if (argc >= 5) {
        L = stod(argv[1]);
        intervals = stoi(argv[2]);
        dt = stod(argv[3]);
        tMax = stod(argv[4]);
    } else if (argc > 1) {
        Usage(argv);
        cout << "Using default parameters\n\n";
    }
    QuantumSystem qSystem(L, intervals, mass);

    int N = enterN();
    
    for (int i = 0; i < N; ++i){
        int n;
        double phi, c;
        cout << "Set state number " << i + 1 << "\n";
        cout << "Set n: ";
        cin >> n;
        cout << "Set phi (phase): ";
        cin >> phi;
        cout << "Set c (amplitude): ";
        cin >> c;
        
        qSystem.addState(n, phi, c);
    }
    qSystem.prepareSystem();
    qSystem.printParams();
    qSystem.simulate(dt, tMax);
}