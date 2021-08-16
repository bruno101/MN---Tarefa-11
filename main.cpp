#include <iostream>
#include <vector>
#include <tuple>
#include <math.h>

using namespace std;

tuple <double, vector<double>> potReg(double n, vector<vector<double>> A, vector<double> v0, double epsilon) {

  double lambdaNovo = 0;
  double lambdaVelho;
  vector<double> vNovo(n);
  vector<double> vVelho(n);
  vector<double> xVelho(n);
  double erro = epsilon+1;

  for (int i = 0; i < n; i++) {
    vNovo[i] = v0[i];
  }

  while (erro > epsilon) {

    lambdaVelho = lambdaNovo;
    double size = 0;
    for (int i = 0; i < n; i++) {
      vVelho[i] = vNovo[i];
      size += vVelho[i]*vVelho[i];
    }
    size = pow(size,0.5);
    for (int i = 0; i < n; i++) {
      vVelho[i] /= size;
    }
    for (int i = 0; i < n; i++) {
      vNovo[i] = 0;
      for (int j = 0; j < n; j++) {
        vNovo[i] += A[i][j]*vVelho[j];
      }
    }
    lambdaNovo = 0;
    for (int i = 0; i < n; i++) {
      lambdaNovo += vVelho[i]*vNovo[i];
    }
    erro = abs((lambdaNovo-lambdaVelho)/lambdaNovo);

  }

  return make_tuple(lambdaNovo, vNovo);

}

vector<vector<double>> calculaInversa (double n, vector<vector<double>> A) {

  vector<vector<double>> X (n, vector<double>(n));

  for (int i = 0; i <= n-1; i++) {
    X[i][i] = 1;
  }

  for (int k = 0; k <= n-1; k++) {

    for (int j = k+1; j <= n-1; j++) {
      X[k][j] = X[k][j]/A[k][k];
      A[k][j] = A[k][j]/A[k][k];
    }
    for (int j = 0; j <= k-1; j++) {
      X[k][j] = X[k][j]/A[k][k];
    }

    X[k][k] = X[k][k]/A[k][k];
    A[k][k] = 1;

    for (int i = 0; i <= n-1; i++) {
      if (i != k) {
        for (int j = k+1; j <= n-1; j++) {
          X[i][j] -= A[i][k]*X[k][j];
          A[i][j] -= A[i][k]*A[k][j];
        }
        for (int j = 0; j <= k-1; j++) {
          X[i][j] -= A[i][k]*X[k][j];
        }
        X[i][k] -= A[i][k]*X[k][k];
        A[i][k] = 0;
      }
    }
  }

  return X;

}

tuple <double, vector<double>> invPower(double n, vector<vector<double>> A, vector<double> v0, double epsilon) {
  vector<vector<double>> invA = calculaInversa(n, A);
  auto res = potReg(n, invA, v0, epsilon);
  return make_tuple(1.0/get<0>(res), get<1>(res));
}

tuple <double, vector<double>> shiftedPower(double n, vector<vector<double>> A, vector<double> v0, double epsilon, double mu) {

  for (int i = 0; i < n; i++) {
    A[i][i] -= mu;
  }
  auto res = invPower(n, A, v0, epsilon);
  return make_tuple(get<0>(res) + mu, get<1>(res));

}



int main() {

  int n;

  cout << "Digite a dimensão da matriz: ";
  cin >> n;

  vector<double> v0(n);
  vector<vector<double>> A(n, vector<double>(n));
  vector<double> v(n);
  double lambda;
  double epsilon;

  for (int i = 0; i < n; i++) {
    v0[i] = 1;
    for (int j = 0; j < n; j++) {
      cout << "Digite o valor de A[" << i+1 << "][" << j+1 << "]: ";
      cin >> A[i][j];
    }
  }

  cout << "Digite o valor da precisão desejada: ";
  cin >> epsilon;

  int e;
  cout << "Escolha o método a ser utilizado:\n1 - Método da potência regular \n2 - Inverse Power Method\n3 - Shifted Power Method\n";
  cin >> e;

  if (e == 1) {

    auto res = potReg(n, A, v0, epsilon);

    cout << "Autovalor dominante: " << get<0>(res) << "\nAutovetor correspondente: (";
    for (int i = 0; i < n; i++) {
      cout << get<1>(res)[i] << ", ";
    }
    cout << "\b\b)";

  } else if (e == 2) {

    auto res = invPower(n, A, v0, epsilon);

    cout << "Autovalor dominante: " << get<0>(res) << "\nAutovetor correspondente: (";
    for (int i = 0; i < n; i++) {
      cout << get<1>(res)[i] << ", ";
    }
    cout << "\b\b)";

  } else {

    double mu;
    cout << "Escolha o número do qual o autovalor deve estar mais próximo:";
    cin >> mu;

    auto res = shiftedPower(n, A, v0, epsilon, mu);

    cout << "Autovalor dominante: " << get<0>(res) << "\nAutovetor correspondente: (";
    for (int i = 0; i < n; i++) {
      cout << get<1>(res)[i] << ", ";
    }
    cout << "\b\b)";

  }

}