#include "GaussInteger.h"
#include<iostream>
GaussInteger::GaussInteger(int Gauss_point_number)
    :Gauss_point_number(Gauss_point_number){

    if (Gauss_point_number == 4) {
        Gauss_coefficient_reference_1D = new std::vector<double>{ 0.3478548451, 0.3478548451, 0.6521451549, 0.6521451549 };
        Gauss_point_reference_1D = new std::vector<double>{ 0.8611363116, -0.8611363116, 0.3399810436, -0.3399810436 };
    }

    else if (Gauss_point_number == 8) {
        Gauss_coefficient_reference_1D = new std::vector<double>{ 0.1012285363, 0.1012285363, 0.2223810345, 0.2223810345, 0.3137066459, 0.3137066459, 0.3626837834, 0.3626837834 };
        Gauss_point_reference_1D = new std::vector<double>{ 0.9602898565, -0.9602898565, 0.7966664774, -0.7966664774, 0.5255324099, -0.5255324099, 0.1834346425, -0.1834346425 };
    }
    else if (Gauss_point_number == 2) {
        Gauss_coefficient_reference_1D = new std::vector < double>{ 1, 1 };
        Gauss_point_reference_1D = new std::vector < double>{ -0.57735026919, 0.57735026919 };
    }
    Gauss_coefficient_local_1D = new std::vector<double>;
    Gauss_point_local_1D = new std::vector<double>;
}
void GaussInteger::generate_Gauss_local_1D(double lowerbound, double upperbound)
{
    Gauss_coefficient_local_1D->clear();
    Gauss_point_local_1D->clear();
    for (int i = 0; i < Gauss_coefficient_reference_1D->size(); i++) {
        Gauss_coefficient_local_1D->push_back((upperbound - lowerbound) * Gauss_coefficient_reference_1D->at(i) / 2);       
    }
    for (int i = 0; i < Gauss_point_reference_1D->size(); i++) {
        Gauss_point_local_1D->push_back((upperbound - lowerbound) * Gauss_point_reference_1D->at(i) / 2 + (upperbound + lowerbound)/2);
    }
}

void GaussInteger::set_bound(double lower, double upper)
{
    lowerbound = lower;
    upperbound = upper;
    generate_Gauss_local_1D(lowerbound,upperbound);
}

double GaussInteger::integral(std::function<double(double)> Func)
{
    generate_Gauss_local_1D(this->lowerbound,this->upperbound);
    
    double result = 0;
    for (int i = 0; i < Gauss_coefficient_local_1D->size(); i++) {
        
        result += Gauss_coefficient_local_1D->at(i) * Func(Gauss_point_local_1D->at(i));
    }
    return result;
}
