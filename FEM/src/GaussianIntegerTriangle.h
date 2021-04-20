#pragma once
#include <vector>
#include <array>
#include <cmath>
#include<functional>
#include <memory>
class GaussianIntegerTriangle
{
public:
    GaussianIntegerTriangle(int Gauss_point_number);

    ~GaussianIntegerTriangle() {
       
    };
    void set_bound(const std::array<std::array<double, 2>, 3>& triangle);
    double integral(std::function<double(std::array<double,2>)> Func);
private:
    int Gauss_point_number;
    std::array<std::array<double,2>,3> triangle_boundary;


    /// Generate the Gauss coefficientsand Gauss points on the reference triangle whose vertices are(0, 0), (1, 0), (0, 1)
    /// Gauss_point_number:the number of Gauss points in the formula.The Gauss formula depends on it.
    /// Gauss_coefficient_reference_triangle, Gauss_point_reference_triangle : the Gauss coefficientsand Gauss points on the reference triangle

    std::unique_ptr< std::vector<double> > Gauss_coefficient_reference_triangle;
    std::unique_ptr< std::vector<std::array<double, 2>> > Gauss_point_reference_triangle;
    std::unique_ptr< std::vector<double>> Gauss_coefficient_local_triangle;
    std::unique_ptr< std::vector<std::array<double, 2>> > Gauss_point_local_triangle;


    //Generate the Gauss coefficients and Gauss points on an arbitrary interval [lower_bound,upper_bound] by using affine tranformation.
    void generate_Gauss_local_triangle(const std::array<std::array<double, 2>, 3>& triangle);

};

