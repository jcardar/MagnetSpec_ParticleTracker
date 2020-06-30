#include <iostream>
#include <string>
#include <fstream>
#include <vector>


void readUnits(std::ifstream &input_stream, std::vector<std::string> &desired_units) {
    std::string length_unit;
    std::string energy_unit;
    std::string angle_unit;
    std::string magnetic_field_unit;

    input_stream >> length_unit >> energy_unit >> angle_unit >> magnetic_field_unit;

    desired_units = {length_unit, energy_unit, angle_unit, magnetic_field_unit};
}

int readNumOf(std::ifstream &input_stream) {
    std::string numStr;
    input_stream >> numStr;
    int numInt = std::stoi(numStr);
    
    return numInt;
}

void readMagnet(std::ifstream &input_stream, int &magNum, std::vector<std::vector<std::vector<double>>> &magInfo) {
    magNum = readNumOf(input_stream);
    std::string tempStr;

    //outside index goes dimensions, position, magnetic field components
    //dimensions go width, length, height
    for(int i=0; i<3; ++i) {
        std::vector<std::vector<double>> tempMagSet;

        for(int j=0; j<magNum; ++j) {
            std::vector<double> tempInfoBits;

            for(int k=0; k<3; ++k) {
                input_stream >> tempStr;
                double tempDbl = std::stod(tempStr);
                tempInfoBits.push_back(tempDbl);
            }
            tempMagSet.push_back(tempInfoBits);
        }
        magInfo.push_back(tempMagSet);
    }
}

void readBeam(std::ifstream &input_stream, std::vector<std::vector<double>> &beamInfo) {
    std::string tempStr;

    //outside index goes position, energy, direction
    for(int i=0; i<3; ++i) {
        std::vector<double> tempInfoBits;

        if(i==1) {
            input_stream >> tempStr;
            double tempDbl = std::stod(tempStr);
            tempInfoBits.push_back(tempDbl);
        }
        else {
            for(int j=0; j<3; ++j) {
                input_stream >> tempStr;
                double tempDbl = std::stod(tempStr);
                tempInfoBits.push_back(tempDbl);
            }
        }
        beamInfo.push_back(tempInfoBits);
    }
}

void readSpread(std::ifstream &input_stream, std::vector<std::vector<double>> &spreadInfo) {
    std::string tempStr;

    //outside index goes position, energy, divergence
    for(int i=0; i<3; ++i) {
        std::vector<double> tempInfoBits;

        if(i==1) {
            input_stream >> tempStr;
            double tempDbl = std::stod(tempStr);
            tempInfoBits.push_back(tempDbl);
        }
        else {
            for(int k=0; k<3; ++k) {
                input_stream >> tempStr;
                double tempDbl = std::stod(tempStr);
                tempInfoBits.push_back(tempDbl);
            }
        }
        spreadInfo.push_back(tempInfoBits);
    }
}

void readScreen(std::ifstream &input_stream, int &screenNum, std::vector<std::vector<std::vector<double>>> &screenInfo) {
    screenNum = readNumOf(input_stream);
    std::string tempStr;

    //outside index goes dimensions, position, angles
    //dimensions go length, height
    //angles go z-axis, x-axis
    for(int i=0; i<3; ++i) {
        std::vector<std::vector<double>> tempScreenSet;

        for(int j=0; j<screenNum; ++j) {
            std::vector<double> tempInfoBits;

            if(i==0) {
                for(int k=0; k<3; ++k) {
                    input_stream >> tempStr;
                    double tempDbl = std::stod(tempStr);
                    tempInfoBits.push_back(tempDbl);
                }
            }
            else {
                for(int k=0; k<2; ++k) {
                    input_stream >> tempStr;
                    double tempDbl = std::stod(tempStr);
                    tempInfoBits.push_back(tempDbl);
                }
            }
            tempScreenSet.push_back(tempInfoBits);
        }
        screenInfo.push_back(tempScreenSet);
    }
}

int main() {
    std::ifstream inputdeck("input_deck.txt");

    if(!inputdeck.is_open()) {
        std::cout << "Input file could not be opened" << std::endl;
        return 1;
    }

    std::vector<std::string> desired_units;
    readUnits(inputdeck, desired_units);

    int number_of_magnets;
    std::vector<std::vector<std::vector<double>>> magnet_information;
    readMagnet(inputdeck, number_of_magnets, magnet_information);

    int number_of_particles = readNumOf(inputdeck);
    std::vector<std::vector<double>> beam_information;
    readBeam(inputdeck, beam_information);

    std::vector<std::vector<double>> beam_spread_information;
    readSpread(inputdeck, beam_spread_information);

    int number_of_screens;
    std::vector<std::vector<std::vector<double>>> screen_information;
    readScreen(inputdeck, number_of_screens, screen_information);

    inputdeck.close();
}