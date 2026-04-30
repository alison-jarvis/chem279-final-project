#include <armadillo>
#include <cmath>

//assuming dipole moment is arma::vec of length 3 = mu
//assuming polarizability tensor is (3 x 3) arma::mat = alpha
double dipole_mag(const arma::vec& mu) {
    return std::sqrt(arma::dot(mu, mu));
}

double iso_polarizability(const arma::mat& alpha){
    return (alpha(0,0) + alpha(1,1) + alpha(2,2)) / 3.0;
}

double aniso_polarizability(const arma::mat& alpha){
    double xy_delta = (alpha(0,0) - alpha(1,1)) * (alpha(0,0) - alpha(1,1));
    double yz_delta = (alpha(1,1) - alpha(2,2)) * (alpha(1,1) - alpha(2,2));
    double xz_delta = (alpha(0,0) - alpha(2,2)) * (alpha(0,0) - alpha(2,2));

    return std::sqrt(0.5 * (xy_delta + yz_delta + xz_delta))
}

int main() {
    return 0;
}
