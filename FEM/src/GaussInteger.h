#pragma once
#include <vector>
#include <cmath>
#include<functional>
class GaussInteger
{
public:
	GaussInteger(int Gauss_point_number);
	
    ~GaussInteger() {
        Gauss_coefficient_reference_1D->clear();
        Gauss_point_reference_1D->clear();
        delete Gauss_coefficient_reference_1D;
        delete Gauss_point_reference_1D;
    };
    void set_bound(double lower, double upper);
    double integral(std::function<double(double)> Func);
private:
    int Gauss_point_number;
    double lowerbound, upperbound;
    /// <summary>
    /// Gauss_coefficient_reference_1D,Gauss_point_reference_1D:the Gauss coefficients and Gauss points on the reference interval [-1,1].
    /// Gauss_coefficient_local_1D,Gauss_point_local_1D: the Gauss coefficients and Gauss points on the local interval.
    /// </summary>
    std::vector<double>* Gauss_coefficient_reference_1D;
    std::vector<double>* Gauss_point_reference_1D;
    std::vector<double>* Gauss_coefficient_local_1D;
    std::vector<double>* Gauss_point_local_1D;

    //Generate the Gauss coefficients and Gauss points on an arbitrary interval [lower_bound,upper_bound] by using affine tranformation.
    void generate_Gauss_local_1D(double lowerbound, double upperbound);
    
};

