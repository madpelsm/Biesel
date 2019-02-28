#include <src/biesel.h>
#include <src/cnoidal.h>
#include <src/stokesII.h>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <map>
#include <string>
std::string filename = "output.txt";

void setParameter(std::string input,
                  std::map<std::string, double>* parameters) {
    std::string param = input.substr(0, input.find("="));
    double argument = std::atof(input.substr(input.find("=") + 1).c_str());
    std::vector<std::string> validCommands = {"H", "T", "d", "x"};
    if (std::find(validCommands.begin(), validCommands.end(), param) !=
        validCommands.end()) {
        parameters->insert(std::pair<std::string, double>(param, argument));
    }
    if (param == "help") {
        printf(
            "Use as <executable> H=1 T=1 d=1 x=15 name=filename.extension\n");
    }
    if (param == "name") {
        filename = input.substr(input.find("=") + 1).c_str();
    }
}

void save_data_to_file(std::string _filename, std::string _data) {
    std::ofstream save_data(_filename);

    if (save_data.is_open()) {
        save_data << _data;
    }
    save_data.close();
}
int main(int argc, char* argv[]) {
    std::map<std::string, double> parameters;
    for (size_t i = 0; i < argc; i++) {
        setParameter(argv[i], &parameters);
    }
    for (auto x : parameters) {
        printf("%s:%f\n", x.first.c_str(), x.second);
    }
    // h H T
    Biesel b1(parameters["d"], parameters["H"], parameters["T"]);
    cnoidal c1(parameters["H"], parameters["T"], parameters["d"]);
    StokesII s1(parameters["H"], parameters["T"], parameters["d"]);

    double x[] = {parameters["x"]};

    std::string data = "";
    int T = 100;
    double dt = 0.1;
    data += "Time [s]\tBiesel [m]\tStokes II \tCnoidal\n";
    for (double t = 0; t < 100; t += dt) {
        // printf("%f \n", (double)b1.eta(x, (double)t));
        data += std::to_string(t) + "\t" + std::to_string(b1.eta(x[0], t)) +
                "\t" + std::to_string(s1.eta(x[0], t)) + "\t" +
                std::to_string(c1.eta(x[0], t)) + "\n";

        // printf("eta(%f,%f) = %f \n", x, (float)t, (double)b1.eta(x,
        // (double)t));
    }

    save_data_to_file(filename, data);
    printf("EOP\n");
    return 0;
}
