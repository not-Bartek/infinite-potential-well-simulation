#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fstream>
#include <algorithm>

#define THROW_ERROR(msg) throw std::runtime_error(std::string(msg) + " [in: " + __FUNCTION__ + " | " + __FILE__ + ":" + std::to_string(__LINE__) + "]")

using namespace std;

enum class SimulationType {
    PotentialWell,
    Oscilator
};

namespace Physics {
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    constexpr long double h = 6.62607015e-34L;
    constexpr long double hbar = h / (2.0L * pi);
}

function<double(double)> makeHermitePolynomial(int n) {
    vector<double> coeffs;
    if (n == 0) {
        coeffs = {1.0};
    } else if (n == 1) {
        coeffs = {0.0, 2.0};
    } else {
        vector<double> prev2 = {1.0};
        vector<double> prev1 = {0.0, 2.0};
        vector<double> curr;

        for(int i = 2; i <= n; ++i) {
            curr.assign(i + 1, 0.0);

            
            for (size_t k = 0; k < prev1.size(); ++k) {
                curr[k + 1] += 2.0 * prev1[k];
            }
            double factor = 2.0 * (i - 1);
            for(size_t k = 0; k < prev2.size(); ++k) {
                curr[k] -= factor * prev2[k];
            }
            prev2 = prev1;
            prev1 = curr;
        }
        coeffs = prev1;
    }

    return [coeffs](double xi) -> double {
        double result = 0.0;
        for (int k = (int)coeffs.size() - 1; k >= 0; --k) {
            result = result * xi + coeffs[k];
        }
        return result;
    };
}

double factorial(int n) {
    double res = 1.0;
    for (int i = 2; i <= n; ++i) res *= i;
    return res;
}

class StationaryState{
    private: 
        const int n;
        const double E;
        vector<double> psi_spatial;

        double CalculateEnergy(SimulationType t, int n, double param, double m) {
            switch(t) {
                case SimulationType::PotentialWell: {
                    double L = param;
                    return (Physics::hbar * n * Physics::pi) * (Physics::hbar * n * Physics::pi) / (2 * m * L * L);
                }
                case SimulationType::Oscilator: {
                    double omega = param;
                    return Physics::hbar * omega * (n + 0.5);
                }
            }
            return 0.0;
        }

        vector<double> CalculatePsiSpatial(SimulationType t, int n, int N, double m, double L = 5.0, double omega = 1.0) {
            vector<double> psi;
            psi.reserve(N);

            switch(t) {
                case SimulationType::PotentialWell: {
                    double dx = L / static_cast<double>(N-1);
                    
                    for(int i = 0; i < N; ++i) {
                        double x = i * dx;
                        psi.push_back(sqrt(2.0/L) * sin(n * Physics::pi / L * x));
                    }
                    break;
                }
                case SimulationType::Oscilator: {
                    
                    double dx = L/ static_cast<double>(N-1);

                    double x_min = -L / 2.0;

                    auto hermit = makeHermitePolynomial(n);
                    double alpha = sqrt(m * omega / Physics::hbar);
                    double norm_coeff = pow(m * omega / (Physics::pi * Physics::hbar), 0.25) / sqrt(pow(2.0, n) * factorial(n));
                    for(int i = 0; i < N; ++i) {
                        double x = x_min + i * dx;
                        double xi = alpha * x;
                        double hval = hermit(xi);
                        psi.push_back(norm_coeff * hval * exp(-xi * xi / 2.0));
                    }
                    break;
                }
            }
            return psi;
        }
    public:

        StationaryState(SimulationType t, int n, int N, double m, double L=0.0, double omega = 0.0) :
            n(n),
            E(t == SimulationType::Oscilator ? CalculateEnergy(t, n, omega, m) : CalculateEnergy(t, n, L, m)),
            psi_spatial(CalculatePsiSpatial(t, n, N, m, L, omega)) {}
        
        

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
        
        SimulationType type;
        const int intervals;
        const double mass;
        const double L;
        const double omega;
    public:
        QuantumSystem(SimulationType type,int intervals, double mass, double L = 0.0, double omega = 0.0):
             type(type), intervals(intervals), mass(mass),  L(L), omega(omega) {}

        void addState(int n, double phi, double c) {
            states.emplace_back(type, n, intervals, mass, L, omega);
            params.emplace_back(phi, c);
        }
        void prepareSystem() {
            StateParams::normalizeVector(params);
        }
        
        int simulate(double dt, double t_max, const string& outputFilename = "simulation.csv") {
            cout << "Starting the simulation...\n";
            ofstream file(outputFilename);
            if (!file.is_open()) return 1;

            vector<complex<double>> psiTotal(intervals);
            int step = 0;
            if(dt <= 0) {
                return 1;
            }
            for (double t = 0.0; t <= t_max; t += dt) {
                fill(psiTotal.begin(), psiTotal.end(), complex<double>(0.0, 0.0));
                
                for(size_t k = 0; k < states.size(); ++k) {
                    const auto& psi_spatial = states[k].getPsi();
                    double E = states[k].getEnergy();
                    double c = params[k].c;
                    double phi = params[k].phi;

                    double omega_t = E / Physics::hbar * t;
                    complex<double> timePhase = polar(c, -(omega_t - phi));

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
            cout << "Simulation finished successfully - saved to " << outputFilename << "\n";
            return 0;
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
    if(N > 0) {
        return N;
    } else {
        return 0;
    }
    
}


void Usage(char* argv[]) {
    cout << "Usage: " << argv[0] << " [L] [N] [dt] [tMax]\n";
}


int main(int argc, char* argv[]) {
    int intervals = 500;
    double dt = 4.0e-7;
    double tMax = 0.005;
    
    double L = 2.0e-6;
    double mass = 4.0e-26;
    double omega = 2.0e4;
    bool animation = 0;
    SimulationType t = SimulationType::Oscilator;
    string outputFilename = "simulation.csv";

    if (argc >= 2) {
        string configFilename = argv[1];
        if (argc >= 3) {
            outputFilename = argv[2];
        }

        ifstream file(configFilename);
        if (!file.is_open()) {
            cerr << "Cannot open config file: " << configFilename << "\n";
            return 1;
        }

        int typeInt, N;
        auto readParam = [&](auto& val, const string& desc) {
            if (!(file >> val)) {
                cerr << "Blod wczytywania parametru: " << desc << "\n";
                exit(1);
            }
            string dummy;
            getline(file, dummy); // zjada resztę linii (komentarze)
        };

        readParam(typeInt, "SimulationType");
        readParam(L, "L");
        readParam(intervals, "intervals");
        readParam(dt, "dt");
        readParam(mass, "mass");
        readParam(tMax, "tMax");
        readParam(omega, "omega");
        readParam(animation, "animation");
        readParam(N, "N (number of states)");
        t = (typeInt == 0) ? SimulationType::PotentialWell : SimulationType::Oscilator;
        QuantumSystem qSystem(t, intervals, mass, L, omega);

        if (N <= 0 || L <= 0 || intervals <= 0 || dt <= 0 || tMax <= 0) {
            cerr << "Incorrect parameters in config file\n";
            return 1;
        }
        for (int i = 0; i < N; ++i) {
            int n;
            double phi, c;
            if (!(file >> n >> phi >> c)) {
                cerr << "Blod wczytywania stanu nr " << i + 1 << "\n";
                return 1;
            }
            string dummy;
            getline(file, dummy);

            if (n + static_cast<int>(t) <= 0 ) {
                cerr << "Incorrect n\n";
                return 1;
            }
            qSystem.addState(n, phi, c);
        }
        file.close();

        qSystem.prepareSystem();
        qSystem.printParams();
        if (qSystem.simulate(dt, tMax, outputFilename)) {
            THROW_ERROR("Symulacja nie powiodla sie");
        }
    } else {
        QuantumSystem qSystem(t, intervals, mass, L, omega);

        int N = enterN();
        if(N <= 0 || L <= 0 || intervals <= 0 || dt <= 0 || tMax <=0 || omega <= 0) {
            return 1;
        }
        for (int i = 0; i < N; ++i){
            int n;
            double phi, c;
            
            cout << "Set state number " << i + 1 << "\n";
            cout << "Set n: ";
            cin >> n;
            if(n < 0) {
                return 1;
            }

            cout << "Set phi (phase): ";
            cin >> phi;
            cout << "Set c (amplitude): ";
            cin >> c;
            
            qSystem.addState(n, phi, c);
        }
        qSystem.prepareSystem();
        qSystem.printParams();
        if(qSystem.simulate(dt, tMax)) {
            THROW_ERROR("Simulation failed");
        }
    }
    if(animation) {
        double x_min = (t == SimulationType::PotentialWell) ? 0 : -L/2.0;
        string command = "python3 visualize.py " + outputFilename + " " + to_string(L) + " " + to_string(x_min);
        int result = std::system(command.c_str());

    if (result != 0) {
        std::cerr << "Animation error\n";
    }
    }
    
    
}