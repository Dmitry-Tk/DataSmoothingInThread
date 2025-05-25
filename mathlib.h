/*! \file MathProcessing.h*/
#ifndef MATHLIB_H
#define MATHLIB_H

#include <iostream>
#include <vector>
#include <QPointF>

/** @defgroup LinearAlgebraFunGroup ФУНКЦИИ ЛИНЕЙНОЙ АЛГЕБРЫ
 *  @{
 */
/// \brief вычисление детерминанта матрицы A размерности \f$N\times N\f$
double getDeterminant(double** A, const uint N);
// Примечание: в случае возникновения сингулярных матриц возвращается false,
// а возвращаемые вектора/матрицы не создаются;
// при успешных вычислениях возвращается true, а возвращаемые вектора/матрицы
// создаются и вычисляются;
// реализованы чистые алгоритмы - входные массивы в аргументах следует создавать перед
// вызовами соответствующих функций
/// \brief решение СЛАУ методом Гаусса
bool Gauss(double* x, double** A, double* b, const uint N);
/// \brief решение СЛАУ методом Гаусса-Жордана
bool Gauss_Jordan(double* x, double** A, double* b, const uint N);
/// \brief решение СЛАУ методом обратной матрицы
bool SLE_InverseMatrix(double* x, double** A, double* b, const uint N);
/// \brief вычисление обратной матрицы
bool getInverseMatrix(double** iA, double** A, const uint N, double& det);
/// \brief решение СЛАУ с симметричной положительно-определенной матрицей
/// с помощью факторизации Холецкого
bool SLE_Cholesky(double* x, double** A, double* b, const uint N);
/// \brief матрично-векторное умножение
void getMatrixVectorProduct(double* result, double** M, double* v, const uint N);
/// \brief матричное умножение
void getMatrixProduct(double** C, double** A, double** B, const uint N);
/// \brief вычисление нижней треугольной матрицы из симметричной
/// положительно-определенной матрицы A через разложение Холецкого A=L*L';
bool Cholesky_Decomposition(double** L, double** A, const uint N);
/// \brief генератор симметричной положительно-определенной
/// матрицы для проверки алгоритма факторизации Холецкого
//void generateSymmPosDefiniteMatrix(double** L, double** A, const uint N);
void generateSymmPosDefiniteMatrix(double** A, const uint N);
/// \brief генератор случайной целочисленной матрицы
void generateRandMatrix(double** M, const uint N);
/// \brief генератор случайной симметричной целочисленной матрицы
void generateSymmMatrix(double** M, const uint N);
/// \brief генератор случайного целочисленного вектора
void generateRandVec(double* v, const uint N);
/// \brief отображение СЛАУ в консоли
void showSLE(double** A, double* b, const uint N);
/// \brief отображение матрицы в консоли
void showMatrix(double** A, const uint N);
/// \warning Возвращаемые объекты (одно/двухмерные массивы находятся на первом позиции в списке аргументов)
/** @} */ // end of LinearAlgebraFunGroup

/** @defgroup InterpolatingFunGroup ИНТЕРПОЛЯЦИОННЫЕ ФУНКЦИИИ
 *  @{
 */
/// \brief интерполяция Лагранжа в точке xi по N точкам (*x,*y)
double interp(double xi, const double* x, const double* y, const uint N);
double LagrangePoly(double xi, const std::vector<double>& x, const std::vector<double>& y);
/// \brief векторная интерполяция в точках *xi по (N+1)-мерному интерполирующему вектору *y(*x)
void interp(double* xi, double* yi, const uint Ni, const uint Nl, const double* x, const double* y, const uint N);
// векторная интерполяция со сглаживанием по N-мерному интерполирующему вектору *y(*x)
// xi - Ni-мерный вектор координат, в которых интерполируются данные
// альтернативная модель интерполяции: по соседям удаленным на расстояния ~+-il*dxi
// ~o+o.......o+o.......o+o.......o+o~ например, здесь по 4 точкам
// . - точки ориг. данных, o - лок.окружение, + - точка интерполяции-среднее от прав-лев соседа, ~ - разрыв
// Условия применения: удобна, если число ориг. данных много (N>Ni),
// а кривую зависимости необходимо представить плавной и сглаженной,
// но малым (меньшим) числом точек - фактически означает сжатие данных
void InterpCompression(double* xi, double* yi, const uint Ni, const uint Nl, const double* x, const double* y, const uint N);
/** @} */ // end of InterpolatingFunGroup

/** @defgroup SmoothingFunGroup ФУНКЦИИ СГЛАЖИВАНИЯ ДАННЫХ
 *  @{
 */
/// \brief полиномиальное сглаживание данных {x, y} минимизацией стандартного отклонения
double PolynomSmoothing(double* xs, double* ys, const uint Ns, double* x, double* y, const uint N, const int n);
// возвращает стандартное отклонение для N-мерного вектора данных,
// сглаженных с помощью полинома оптимального порядка, но не выше порядка n+1,
// а также возвращает коэффициенты этого сглаживающего полинома a[0]+a[1]*x+a[2]*x^2+...a[n]*x^n
// и его порядок;
double getCoeffOptimalSmoothPoly(std::vector<double>& a, uint n, const double* x, const double* y, const uint N);
// для случая полинома фиксированного порядка (a.size()-1)
double getCoeffSmoothPoly(std::vector<double>& a, const double* x, const double* y, const uint N);
double getCoeffSmoothPoly(std::vector<double>& a, const std::vector<QPointF> points);
// расчёт матрицы X и вектора f, используемых для нахождения коэффициентов сглаживающего полинома
void getCoeffMatrix(double** X, double* f, const uint n, const double* x, const double* y, const uint N);
void getCoeffMatrix(double** X, double* f, const uint n, const std::vector<QPointF>& data);
/// \brief возвращает значение полинома n-го порядка с вектором коэффициентов a
double Polynomial(double x, const std::vector<double>& a);
/// \brief возвращает значение полинома n-го порядка с вектором корней roots
double PolynomialR(double x, const std::vector<double>& roots);
/** @} */ // end of SmoothingFunGroup

/** @defgroup ReadingWritingFileFunGroup ФУНКЦИИ ЗАПИСИ/ЧТЕНИЯ В/ИЗ ФАЙЛА
 * @{
 */
void WriteFile(const char* FileName, double* x, double* y, const int N);
void ReadFile(const char* FileName, double* x, double* y, int& N);
void WriteFile(const char* FileName, const std::vector<double>& x, const std::vector<double>& y);
void ReadFile(const char* FileName, std::vector<double>& x, std::vector<double>& y);
//void WriteFile(const char* FileName, const std::vector<QPointF>& data);
//void ReadFile(const char* FileName, std::vector<QPointF>& data);
/** @}*/ // end of ReadingWritingFileFunGroup

/** @defgroup AuxFunGroup ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
 *  @{
 */
// генератор модельных данных
double* generateData(const uint N, double xl, double xr, double *x);
// генератор случайных чисел rnd() = [-1.0...1.0]
double rnd();
// генератор случайной линейной функции
double RndFunction(const double xi, const bool boolJumps);
// измеритель времени в тактах процессора
//static inline unsigned long long int tick();
/** @} */ // end of AuxFunGroup

/** @defgroup NodNokGroup НОД и НОК
 * @{
 */
void Factorization(long long int f, std::vector<int>& Factor);
long long int getNOD(std::vector<int>& Factor1, std::vector<int>& Factor2);
long long int getNOD(long long int f1, long long int f2);

#endif // MATHLIB_H
