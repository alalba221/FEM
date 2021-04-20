#include "GaussianIntegerTriangle.h"

GaussianIntegerTriangle::GaussianIntegerTriangle(int Gauss_point_number)
    :Gauss_point_number(Gauss_point_number)
{

    if (Gauss_point_number == 3) {
        Gauss_coefficient_reference_triangle = std::make_unique<std::vector<double>, std::initializer_list<double>>(
            { 1.0 / 6.0,1.0 / 6.0,1.0 / 6.0 }
        );
        Gauss_point_reference_triangle = std::make_unique<std::vector<std::array<double, 2>>, std::initializer_list<std::array<double, 2>>>(
            { {1.0 / 2.0,0},{1.0 / 2.0,1.0 / 2.0},{0, 1.0 / 2.0} }
        );
    }
    else if (Gauss_point_number == 4) {
        Gauss_coefficient_reference_triangle = std::make_unique<std::vector<double>, std::initializer_list<double>>(
            { 
                (1.0 - 1.0 / sqrt(3.0)) / 8.0, 
                (1.0 - 1.0 / sqrt(3.0)) / 8.0, 
                (1.0 + 1.0 / sqrt(3.0)) / 8.0, 
                (1.0 + 1.0 / sqrt(3.0)) / 8.0 
            }
        );
        Gauss_point_reference_triangle = std::make_unique<std::vector<std::array<double, 2>>, std::initializer_list<std::array<double, 2>>>(
            {
                {(1.0 / sqrt(3.0) + 1.0) / 2.0, (1.0 - 1.0 / sqrt(3.0)) * (1.0 + 1.0 / sqrt(3.0)) / 4.0},
                {(1.0 / sqrt(3.0) + 1.0) / 2.0, (1.0 - 1.0 / sqrt(3.0)) * (1.0 - 1.0 / sqrt(3.0)) / 4.0},
                {(-1.0 / sqrt(3.0) + 1.0) / 2.0, (1.0 + 1.0 / sqrt(3.0)) * (1.0 + 1.0 / sqrt(3.0)) / 4.0},
                {(-1.0 / sqrt(3.0) + 1.0) / 2.0, (1.0 + 1.0 / sqrt(3.0)) * (1.0 - 1.0 / sqrt(3.0)) / 4.0 }
            }
        );
    }
    else if (Gauss_point_number == 9) {
        Gauss_coefficient_reference_triangle = std::make_unique<std::vector<double>, std::initializer_list<double>>(
            {
                64.0 / 81.0 * (1.0 - 0.0) / 8.0, 
                100.0 / 324.0 * (1.0 - sqrt(3.0 / 5.0)) / 8.0, 
                100.0 / 324.0 * (1.0 - sqrt(3.0 / 5.0)) / 8.0, 
                100.0 / 324.0 * (1.0 + sqrt(3.0 / 5.0)) / 8.0, 
                100.0 / 324.0 * (1.0 + sqrt(3.0 / 5.0)) / 8.0, 
                40.0 / 81.0 * (1.0 - 0.0) / 8.0, 
                40.0 / 81.0 * (1.0 - 0.0) / 8.0, 
                40.0 / 81.0 * (1.0 - sqrt(3.0 / 5.0)) / 8.0, 
                40.0 / 81.0 * (1.0 + sqrt(3.0 / 5.0)) / 8.0
            }
        );
        Gauss_point_reference_triangle = std::make_unique<std::vector<std::array<double, 2>>, std::initializer_list<std::array<double, 2>>>(
            {
                {(1.0 + 0.0) / 2.0, (1.0 - 0.0) * (1.0 + 0.0) / 4.0},
                {(1.0 + sqrt(3.0 / 5.0)) / 2.0, (1.0 - sqrt(3.0 / 5.0)) * (1.0 + sqrt(3.0 / 5.0)) / 4.0},
                {(1.0 + sqrt(3.0 / 5.0)) / 2.0, (1.0 - sqrt(3.0 / 5.0)) * (1.0 - sqrt(3.0 / 5.0)) / 4.0},
                {(1.0 - sqrt(3.0 / 5.0)) / 2.0, (1.0 + sqrt(3.0 / 5.0)) * (1.0 + sqrt(3.0 / 5.0)) / 4.0},
                {(1.0 - sqrt(3.0 / 5.0)) / 2.0, (1.0 + sqrt(3.0 / 5.0)) * (1.0 - sqrt(3.0 / 5.0)) / 4.0},
                {(1.0 + 0.0) / 2.0, (1.0 - 0.0) * (1.0 + sqrt(3.0 / 5.0)) / 4.0},
                {(1.0 + 0.0) / 2.0, (1.0 - 0.0) * (1.0 - sqrt(3.0 / 5.0)) / 4.0},
                {(1.0 + sqrt(3.0 / 5.0)) / 2.0, (1.0 - sqrt(3.0 / 5.0)) * (1.0 + 0.0) / 4.0},
                {(1.0 - sqrt(3.0 / 5.0)) / 2.0, (1.0 + sqrt(3.0 / 5.0)) * (1.0 + 0.0) / 4.0}
            }
        );
    }

    Gauss_coefficient_local_triangle = std::make_unique < std::vector<double>>();
    Gauss_point_local_triangle = std::make_unique< std::vector< std::array<double, 2> > >();
}

void GaussianIntegerTriangle::set_bound(const std::array<std::array<double, 2>, 3>& triangle)
{
    triangle_boundary = triangle;
    generate_Gauss_local_triangle(this->triangle_boundary);
}

double GaussianIntegerTriangle::integral(std::function<double(std::array<double, 2>)> Func)
{
    double result = 0;
    for (int i = 0; i < Gauss_coefficient_local_triangle->size(); i++) {

        result += Gauss_coefficient_local_triangle->at(i) * Func(Gauss_point_local_triangle->at(i));
    }
    return result;
}

void GaussianIntegerTriangle::generate_Gauss_local_triangle(const std::array<std::array<double, 2>, 3>& triangle)
{
    Gauss_coefficient_local_triangle->clear();
    Gauss_point_local_triangle->clear();
    
    double x0 = triangle[0][0];
    double y0 = triangle[0][1];
    double x1 = triangle[1][0];
    double y1 = triangle[1][1];
    double x2 = triangle[2][0];
    double y2 = triangle[2][1];
    
    double Jacobi = fabs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));

    for (int i = 0; i < Gauss_coefficient_reference_triangle->size(); i++) {
        Gauss_coefficient_local_triangle->push_back(Gauss_coefficient_reference_triangle->at(i) * Jacobi);
    }

    for (int i = 0; i < Gauss_point_reference_triangle->size(); i++) {
        
        std::array<double, 2> point;

        point[0] = x0 + (x1 - x0) * Gauss_point_reference_triangle->at(i)[0] + (x2 - x0) * Gauss_point_reference_triangle->at(i)[1];
        point[1] = y0 + (y1 - y0) * Gauss_point_reference_triangle->at(i)[0] + (y2 - y0) * Gauss_point_reference_triangle->at(i)[1];
        Gauss_point_local_triangle->push_back(point);
    }
}
