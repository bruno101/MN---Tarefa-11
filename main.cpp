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
      xVelho[i] = vVelho[i]/size;
    }
    for (int i = 0; i < n; i++) {
      vNovo[i] = 0;
      for (int j = 0; j < n; j++) {
        vNovo[i] += A[i][j]*xVelho[j];
      }
    }
    lambdaNovo = 0;
    for (int i = 0; i < n; i++) {
      lambdaNovo += xVelho[i]*vNovo[i];
    }

    erro = abs((lambdaNovo-lambdaVelho)/lambdaNovo);
    if (lambdaNovo == 0) {
      if (abs(lambdaNovo-lambdaVelho) < erro) {
        erro = 0;
      } else {
        erro = epsilon+1;
      }
    }

  }

  return make_tuple(lambdaNovo, xVelho);

}

tuple <vector<vector<double>>, vector<vector<double>>> decompLU (double n, vector<vector<double>> A) {

  vector<vector<double>> L(n, vector<double>(n)), U(n, vector<double>(n));
  double S;

  for (int j = 0; j < n; j++) {

    for (int i = 0; i <= j; i++) {
      U[i][j] = A[i][j];
      for (int k = 0; k < i; k++) {
        U[i][j] -= L[i][k]*U[k][j];
      }
    }

    L[j][j] = 1;

    for (int i = j+1; i < n; i++) {
      L[i][j] = A[i][j];
      for (int k = 0; k < j; k++) {
        L[i][j] -= L[i][k]*U[k][j];
      }
      L[i][j] /= U[j][j];
    }

  }

  return make_tuple(L, U);

}

vector<double> LUsolver (double n, vector<vector<double>> L, vector<vector<double>> U, vector<double> b) {

  vector<double> x(n), y(n);

  for (int i = 0; i < n; i++) {
    y[i] = b[i];
    for (int k = 0; k < i; k++) {
      y[i] -= L[i][k]*y[k];
    }
  }

  for (int i = n-1; i >= 0 ; i--) {
    x[i] = y[i];
    for (int k = i+1; k < n; k++) {
      x[i] -= U[i][k]*x[k];
    }
    x[i] /= U[i][i];
  }

  return x;

}


tuple <double, vector<double>> invPower(double n, vector<vector<double>> A, vector<double> v0, double epsilon) {

  auto LU = decompLU(n, A);

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
      xVelho[i] = vVelho[i]/size;
    }
    vNovo = LUsolver(n, get<0>(LU), get<1>(LU), xVelho);
    lambdaNovo = 0;
    for (int i = 0; i < n; i++) {
      lambdaNovo += xVelho[i]*vNovo[i];
    }

    erro = abs((lambdaNovo-lambdaVelho)/lambdaNovo);
    if (lambdaNovo == 0) {
      if (abs(lambdaNovo-lambdaVelho) < erro) {
        erro = 0;
      } else {
        erro = epsilon+1;
      }
    }

  }

  return make_tuple(1.0/lambdaNovo, xVelho);

}

tuple <double, vector<double>> shiftedPower(double n, vector<vector<double>> A, vector<double> v0, double epsilon, double mu) {

  for (int i = 0; i < n; i++) {
    A[i][i] -= mu;
  }
  auto res = invPower(n, A, v0, epsilon);
  return make_tuple(get<0>(res) + mu, get<1>(res));

}

tuple<vector<double>, vector<vector<double>>> shifterPowerMultipleValues(double n, vector<vector<double>> A, vector<double> v0, double epsilon, int s) {

  vector<double> autovalores;
  vector<vector<double>> autovetores;

  auto menorSol = invPower(n, A, v0, epsilon);
  double menorAutovalor = get<0>(menorSol);
  auto maiorSol = potReg(n, A, v0, epsilon);
  double maiorAutovalor = get<0>(maiorSol);

  if (maiorAutovalor != maiorAutovalor) {
    cout << "Potência regular retornou nan.\n";
    return make_tuple(autovalores, autovetores);
  } else if (menorAutovalor != menorAutovalor) {
    cout << "Potência inversa retornou nan.\n";
    return make_tuple(autovalores, autovetores);
  }

  if (abs(maiorAutovalor - menorAutovalor) < 0.0001) {
    autovalores.push_back(get<0>(menorSol));
    autovetores.push_back(get<1>(menorSol));
    return make_tuple(autovalores, autovetores);
  }
  if (abs(maiorAutovalor + menorAutovalor) < 0.0001) {
    autovalores.push_back(get<0>(menorSol));
    autovetores.push_back(get<1>(menorSol));
    autovalores.push_back(get<0>(maiorSol));
    autovetores.push_back(get<1>(maiorSol));
    return make_tuple(autovalores, autovetores);
  }

  double intervalo = (abs(maiorAutovalor)-abs(menorAutovalor))/s;
  //autovalores negativos
  double mu = -abs(maiorAutovalor)-intervalo;
  double autovalorAnterior = -abs(maiorAutovalor)-1000;
  while (1) {

    mu += intervalo;
    auto sol = shiftedPower(n, A, v0, epsilon, mu);
    if (get<0>(sol) - autovalorAnterior > 0.0001) {
      autovalores.push_back(get<0>(sol));
      autovetores.push_back(get<1>(sol));
      autovalorAnterior = get<0>(sol);
    }

    if (mu > -abs(menorAutovalor)) {
      break;
    }

  }
  //autovalores positivos
  mu = max(abs(menorAutovalor)-intervalo, mu);
  while (1) {

    mu += intervalo;
    auto sol = shiftedPower(n, A, v0, epsilon, mu);
    if (get<0>(sol) - autovalorAnterior > 0.0001) {
      autovalores.push_back(get<0>(sol));
      autovetores.push_back(get<1>(sol));
      autovalorAnterior = get<0>(sol);
    }

    if (mu > abs(maiorAutovalor)) {
      break;
    }

  }

  return make_tuple(autovalores, autovetores);

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

  while (1) {

  int e;
  cout << "Escolha o método a ser utilizado:\n1 - Método da potência regular \n2 - Inverse Power Method\n3 - Shifted Power Method (encontrar autovalor próximo a um valor específico)\n4 - Shifted Power Method (encontrar múltiplos autovalores)\n";
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

    cout << "Autovalor: " << get<0>(res) << "\nAutovetor correspondente: (";
    for (int i = 0; i < n; i++) {
      cout << get<1>(res)[i] << ", ";
    }
    cout << "\b\b)";

  } else if (e == 3) {

    double mu;
    cout << "Escolha o número do qual o autovalor deve estar mais próximo:";
    cin >> mu;

    auto res = shiftedPower(n, A, v0, epsilon, mu);

    cout << "Autovalor: " << get<0>(res) << "\nAutovetor correspondente: (";
    for (int i = 0; i < n; i++) {
      cout << get<1>(res)[i] << ", ";
    }
    cout << "\b\b)";

  } else {

    int s;
    cout << "Escolha o número de partes em que se deve dividir os intervalos: ";
    cin >> s;
    cout << "\n";
    auto res = shifterPowerMultipleValues(n, A, v0, epsilon, s);

    vector<double> autovalores = get<0>(res);
    vector<vector<double>> autovetores = get<1>(res);
    cout << "Pares de autovalores e autovetores encontrados:\n";
    for (int i = 0; i < size(autovalores); i++) {
      cout << autovalores[i] << ", (";
      for (int j = 0; j < n; j++) {
        cout << autovetores[i][j] << ", ";
      }
      cout << "\b\b)\n";
    }

  }

  int e1;
  cout << "\n\nVocê deseja continuar?\n1 - Sim\n2 - Não\n";
  cin >> e1;
  if (e1 == 2) {
    break;
  }
  cout << "\n";

  }

}