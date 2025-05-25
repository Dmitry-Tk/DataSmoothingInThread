/*! \file MathProcessing.cpp*/

#include "mathlib.h"
#include <iostream>
#include <fstream>
#include <cstdlib>  // random()
#include <cmath>    // exp(), M_PI
#include <stdio.h>  // FILE, fopen(), fclose(), fprintf()
#define SLINE "--------"

// возникла неоднозначность функции isnan т.к. она определена в двух местах
// using namespace std;

double Amp = 10.0;// амплитуда случайных чисел в генераторах случайных матриц/векторов
extern int N;// длина векторов данных
//const short int Nl = N/10;// длина секций, на которых параметры RndFunction сохраняются
const short int Nl = 100/10;// длина секций, на которых параметры RndFunction сохраняются
extern const double yTop = 3;
//extern const double yBottom = -3;

///
/// \details Точечная интерполяция (в точке xi) по N точкам {*x,*y}.
/// Основана на интерполяционном многочлене Лагранжа N-го порядка.
/// \param[in] xi - координаты точек интерполяции
/// \param[in] {x, y} - вектора данных по которым выполняется интерполяция
/// \return интерполированное значение
///
double interp(double xi, const double* x, const double* y, const uint N)
{
    double yi = 0.0, P;
    for (uint n = 0; n < N; ++n)
    {
        P = 1.0;
        for (uint k = 0; k < N; ++k)
            if (k != n)
                P *= (xi-x[k])/(x[n]-x[k]);
        yi += y[n]*P;
    }

    return yi;
}
///
/// \brief LagrangePoly - полином Лагранжа
/// \param xi
/// \param x
/// \param y
/// \return
///
double LagrangePoly(double xi, const std::vector<double>& x, const std::vector<double>& y)
{
    double S = 0.0, P;
    for (size_t n = 0; n < x.size(); n++)
    {
        P = 1.0;
        for (size_t k = 0; k < x.size(); k++)
            if (k != n)
                P *= (xi - x[k])/(x[n] - x[k]);
        S += y[n]*P;
    }
    return S;
}
///
/// \details Векторная интерполяция по N-мерному интерполирующему вектору *y(*x)
/// xi - Ni-мерный вектор координат, в которых интерполируются данные.
/// Модель интерполяции - по ближайшим соседям:
///     \f$\sim....oo+oo...\sim...oo+oo...\sim...oo+oo...\sim\f$ например, здесь по 4 точкам
/// \f$.\f$ - точки ориг. данных, \f$o\f$ - лок.окружение, \f$+\f$ - точка интерполяции, \f$\sim\f$ - разрыв
/// \param[out] xi - вектор абсцисс интерполированных данных
/// \param[out] yi - вектор интерполированных данных
/// \param[in] Ni - число интерполирующих точек
/// \param[in] Nl - число точек используемых для интерполяции в одной точке
/// \param[in] x  - вектор абсцисс исходных данных *y(*x) подразумевается (N+1)-мерным [0..N]
/// \param[in] y  - вектор исходных данных *y(*x) подразумевается (N+1)-мерным [0..N]
/// \param[in] N  - размерность исходных данных
/// \warning Условия применения: удобна, если число оригинальных данных мало (N<Ni или даже N<<Ni),
/// а контур кривой зависимости необходим более плавный.
///
void interp(double* xi, double* yi, const uint Ni, const uint Nl, const double* x, const double* y, const uint N)
{
    double* xx = new double [Nl];// локальные данные, используемые при интерполяции
    double* yy = new double [Nl];// локальные данные, используемые при интерполяции
    double  dxi = (x[N-1] - x[0])/(Ni-1);// шаг интерполяции
    int index = 0;// индекс ближайшей точки к интерполируемой

    xi[0] = x[0];       yi[0] = y[0];// левая крайняя точка
    xi[Ni-1] = x[N-1];  yi[Ni-1] = y[N-1];// правая крайняя точка
    for (uint j = 1; j < Ni-1; ++j)
    {
        xi[j] = x[0] + j*dxi;// координаты равноотстоящих интерполируемых точек
        // модель интерполяции: по ближайшим соседям
        // [j-1][j]i[j+1][j+1]
        // ~....oo+oo...~
        while (xi[j] > x[index] && (index <= static_cast<int>(N)-1))
            index++;
        for (uint ii = 0; ii < Nl/2; ++ii)
            --index;
        if (index < 0) index = 0;
        if (index > static_cast<int>(N)-1) index = (static_cast<int>(N)-1) - (static_cast<int>(Nl)-1);
        for (uint il= 0; il < Nl; ++il)
        {
            xx[il] = x[static_cast<uint>(index)+il];
            yy[il] = y[static_cast<uint>(index)+il];
        }
        yi[j] = interp(xi[j], xx, yy, Nl);
    }
    delete[] xx;
    delete[] yy;
}
///
/// \details Векторная интерполяция со сглаживанием по N-мерному интерполирующему вектору *y(*x)
/// xi - Ni-мерный вектор координат, в которых интерполируются данные.
/// Альтернативная модель интерполяции: по соседям удаленным на расстояния \f$\sim+-il*dxi\f$
/// \f$\sim o+o.......o+o.......o+o.......o+o\sim\f$ например, здесь по 4 точкам
/// Обозначения: \f$.\f$ - точки ориг. данных; \f$o\f$ - лок.окружение;
/// \f$+\f$ - точка интерполяции-среднее от прав-лев соседа; \f$\sim\f$ - разрыв.
/// \param[out] xi - вектор абсцисс интерполированных данных
/// \param[out] yi - вектор интерполированных данных
/// \param[in] Ni - число интерполирующих точек
/// \param[in] Nl - число точек используемых для интерполяции в одной точке
/// \param[in] x  - вектор абсцисс исходных данных *y(*x) подразумевается (N+1)-мерным [0..N]
/// \param[in] y  - вектор исходных данных *y(*x) подразумевается (N+1)-мерным [0..N]
/// \param[in] N  - размерность исходных данных
/// \warning Условия применения: удобна, если число ориг. данных много (N>Ni), а кривую зависимости
/// необходимо представить плавной и сглаженной, но малым (меньшим) числом точек - фактически
/// означает сжатие данных.
void InterpCompression(double* xi, double* yi, const uint Ni, const uint Nl, const double* x, const double* y, const uint N)
{
    // yi,xi вычисляемые вектора интерполированных данных
    // Nl - число точек используемое для интерполяции в одной точке
    double* xx = new double [Nl];// локальные данные, используемые при интерполяции
    double* yy = new double [Nl];// локальные данные, используемые при интерполяции
    double  dxi = (x[N-1] - x[0])/(Ni-1);// шаг интерполяции
    int index = 0, index0 = 0;// индекс ближайшей к интерполируемой (справа) точки
    double Dx = static_cast<double>(Nl-1)/2*dxi;//полуширина интерполяционного интервала
    double xc;// центральная точка - точка интерполяции

    yi[0] = y[0]; yi[Ni-1] = y[N-1];// левая и правая крайние точки
    // генерируем вектор координат равноотстоящих интерполируемых точек
    for (uint j = 0; j < Ni; ++j)
        xi[j] = x[0] + j*dxi;
    for (uint j = 1; j < Ni-1; ++j)
    {
        // альтернативная модель интерполяции: по ближайшим соседям,
        // но равноудаленным на расстояния ~+-il*dxi
        // ~i1[j]i1+1~  ~i2[j+1]i2+2~  ~i3[j+2]i3+1  ~i4[j+3]i4+1~
        // ~..o+o...........o+o............o+o...........o+o~
        // ~..o+o...........o+o...........o+o............o+o...........o+o~
        xc = xi[j];//центральная точка - точка интерполяции
        // смещаем центры интерполяционных интервалов на краях всей области
        // интерполяции во внутрь, чтобы не выйти за границы диапазона координаты x
        if (xc-Dx < xi[0])
            xc = xi[0] + 1.01*Dx;
        if (xc+Dx > xi[Ni-1])
            xc = xi[Ni-1] - 1.01*Dx;
        // cout <<"j = "<<j<<" : "<<xi[0]<<"\t"<<xi[Ni-1]<< endl;
        // генерируем локальное окружение для текущей точки интерполяции
        for (uint il = 0; il < Nl; ++il)
        {
            while (xc-Dx+il*dxi > x[index])
                index++;
            if (il == 0)
                index0 = index;
            xx[il] = 0.5*(x[index] + x[index-1]);
            yy[il] = 0.5*(y[index] + y[index-1]);
            //cout <<"\t"<<xx[il]<<"\t"<<yy[il]<< endl;
        }
        yi[j] = interp(xi[j], xx, yy, Nl);
        index = index0;
    }
    delete[] xx;
    delete[] yy;
}
//
// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
//
short int _variant = 0;// выбор модельной функции 0...4 в генераторе данных
///
/// \brief Генератор данных
///
double* generateData(const uint N, double xl, double xr, double *x)
{
    double* y = new double [N+1];
    double dx = (xr - xl)/N;

    for (uint i = 0; i < N+1; ++i)
    {
        x[i] = xl + i*dx;
        switch (_variant)
        {
        case 0: y[i] = exp(-0.25*x[i]*x[i])*sin(10*M_PI*x[i]/xr);
            break;
        case 1: y[i] = exp(-0.25*x[i]*x[i])*cos(10*M_PI*x[i]/xr);
            break;
        case 2: y[i] = exp(-1.0*fabs(x[i]))*sin(10*M_PI*x[i]/xr);
            break;
        case 3: y[i] = exp(-1.0*fabs(x[i]))*cos(10*M_PI*x[i]/xr);
            break;
        default:
        {
            if (x[i] != 0.0)
                y[i] = sin(2*M_PI*x[i])/(2*M_PI*x[i]);
            else
                y[i] = 1;
        }
            break;
        }
    }
    return y;
}
///
/// \details Генератор случайных чисел rnd() = [-1.0...1.0].
///
double rnd()
{
    return 2*static_cast<double>(rand())/RAND_MAX - 1.0;
}
///
/// \details Генератор линейной функции, меняющей случайным образом свои параметры
///
double RndFunction(const double xi, const bool boolJumps)
{
    // переменные и параметры управления генерацией
    // псевдослучайной ф-ции
    static int jl = 0;
    static double a = rnd();// наклон
    static double b = rnd();// отсечка
    static double xlc = 0;  // текущее начало координат
    double yi;
    if (xi == 0.0)
    {
        xlc = xi;
        jl = 0;
    }
    if (++jl == Nl)
    {
        jl = 0;
        // модифицируем b, если со скачками
        // плавный вариант - b = yi(xi)
        if (boolJumps)
            b = rnd();
        else
            b = a*(xi - xlc) + b;
        xlc = xi;
        a = rnd();
    }
    yi = a*(xi - xlc) + b;
    if (fabs(yi) > 2*yTop)
    {
        b = 0;
        return b;
    }
    // ограничиваем кривую полосой |yi|<yTop
    if (yi > yTop)
    {
        a *= -1;
        xlc = xi;
        yi = yTop - (yi - yTop);
    }
    else
        if (yi < -yTop)
        {
            a *= -1;
            xlc = xi;
            yi = -yTop - (yTop + yi);
        }
    return yi;
}
/*************************************************************
 * измеритель времени в тактах процессора                    *
 * 2.50 GHz - тактовая частота процессора                    *
 * пример измерений                                          *
 * --------------------------------------------------------- *
 *    double sec = 1.0, ms = sec/1000, mks = ms/1000;//[sec] *
 *    double ut = ms;                                        *
 *    double GHz = 1000000000;//[Hz]                         *
 *    double T = 1/(2.5*GHz);                                *
 *    tstart = tick();                                       *
 *    int s = 1;                                             *
 *    for (int i= 1; i<277000000; i++)//~ 1 sec               *
 *        s *= i;                                            *
 *    tstop = tick();                                        *
 *    if (tstop > tstart) ticks = tstop - tstart;            *
 *    t = ticks*T/ut;//время выполнения в ms                 *
 * --------------------------------------------------------- *
 * если измерять с точностью до ~ms, то можно                *
 * использовать функцию clock() <ctime>                      *
 * с помощью переменных типа clock_t                         *
 *************************************************************/
/*inline unsigned long long int tick()//static
{
    unsigned long long int d;
    __asm__ __volatile__ ("rdtsc" : "=A" (d));
    return d;
}
*/

///
/// \details Вычисление детерминанта матрицы A размерности \f[N\times N\f]
/// \param[in] A - матрица
/// \param[in] N - размерность матрицы
/// \return определитель матрицы A
///
double getDeterminant(double** A, const uint N)
{    // приведение матрицы A к треугольному виду
    double det = 1.0;//детерминант
    double sign = 1.0;
    // подразумеваемый формат матрицы: A[row][column]
    for (uint m = 0; m < N-1; m++)//m номер опорной строки или диаг. элемента
    {
        if (A[m][m] == 0.0)//перестановка строк
        {
            uint mm = m;
            while (mm < N && A[mm][m] == 0.0)
                mm++;
            if (mm == N)//т.е. матрица сингулярная
            {
                std::cout <<"The matrix is a singular ! det = 0"<< std::endl;
                return 0.0;
            }
            double var;
            for (uint i = 0; i < N; i++)
            {
                var = A[m][i];
                A[m][i] = A[mm][i];
                A[mm][i] = var;
            }
            sign *= -1.0;// перестановка любых двух строк/столбцов меняет знак определителя
        }
        det *= A[m][m];
        for (uint r = m+1; r < N; r++)//r номер текущей строки; пробегаем по всем строкам m+1..N-1
        {
            double q = A[r][m]/A[m][m];
            for (uint c = m; c < N; c++)//c номер столбца
                A[r][c] -= q*A[m][c];//A'[r][c] = A[r][c] - q*A[m][c];
        }
    }
    det *= sign*A[N-1][N-1];
    return det;
}
///
/// \details Решение системы линейных алгебраических уравнений (СЛАУ): \f[A\cdot x = b \f]
/// классическим методом Гаусса - приведением к треугольному виду,
/// а затем методом исключения неизвестных.
/// Классический универсальный метод, основное достоинство - скорость.
/// \param[out] x - вектор корней - решение СЛАУ
/// \param[in] A - матрица
/// \param[in] b - вектор свободных членов
/// \param[in] N - размерность СЛАУ
/// \return bool - true/false - успешные/неуспешные вычисления
/// \warning Матрица A и вектор свободных членов b в процессе расчёта изменяются !
///
bool Gauss(double* x, double** A, double* b, const uint N)
{
    double det = 1.0;// детерминант
    // подразумеваемый формат матрицы: A[row][column]
    // приведение матрицы A к треугольному виду
    for (uint m = 0; m < N-1; m++)//m номер опорной строки или диаг. элемента
    {
        if (A[m][m] == 0.0)//перестановка строк
        {
            uint mm = m;
            while (mm < N && A[mm][m] == 0.0)
                mm++;
            if (mm == N)//т.е. матрица сингулярная
            {
                std::cout <<"The matrix is a singular ! det = 0"<<std::endl;
                return false;
            }
            double helpvar;
            for (uint i = 0; i < N; i++)
            {
                helpvar = A[m][i];
                A[m][i] = A[mm][i];
                A[mm][i] = helpvar;
            }
            helpvar = b[m];
            b[m] = b[mm];
            b[mm] = helpvar;
        }
        det *= A[m][m];
        for (uint r = m+1; r < N; r++)//k номер текущей строки; пробегаем по всем строкам m+1..N-1
        {
            double q = A[r][m]/A[m][m];
            for (uint c = m; c < N; c++)//n номер столбца
                A[r][c] -= q*A[m][c];//A'[r][c] = A[r][c] - q*A[m][c];
            b[r] -= q*b[m];//b'[r] = b[r] - q*b[m];
        }
        // отображение СЛАУ после текущей итерации
        // cout <<"m = "<<m<< endl;
        // showSLE(A,b,N);
    }
    det *= A[N-1][N-1];
    //cout <<"det = "<<det<< endl;
    if (det == 0.0 || std::isnan(det))
    {
        std::cout <<"The matrix is a singular!: det = "<<det<<std::endl;
        return false;
    }
    // находим решение последовательно искючая неизвестные
    for (int r = static_cast<int>(N)-1; r > -1; r--)//номер текущей строки
    {
        double S = 0;
        for (int c = static_cast<int>(N)-1; c > r; c--)//номер текущего столбца
            S += A[r][c]*x[c];
        x[r] = (b[r] - S)/A[r][r];
    }
    return true;
}
///
/// \details Решение системы линейных алгебраических уравнений (СЛАУ)
/// методом Гаусса - Жордана (обращения матрицы).
/// \param[out] x - вектор корней - решение СЛАУ
/// \param[in] A - матрица
/// \param[in] b - вектор свободных членов
/// \param[in] N - размерность СЛАУ
/// \return bool - true/false - успешные/неуспешные вычисления
/// \warning Матрица A и вектор свободных членов b в процессе расчёта изменяются !
///
bool Gauss_Jordan(double* x, double** A, double* b, const uint N)
{
    double det = 1.0;// детерминант
    // подразумеваемый формат матрицы: A[row][column]
    // приведение матрицы A к диагональной единичной
    for (uint m = 0; m < N; m++)//m номер опорной строки или диаг. элемента
    {
        //перестановка строк, если опорный элемент нулевой
        if (A[m][m] == 0.0)
        {
            uint mm = m;
            while (mm < N && A[mm][m] == 0.0)
                mm++;
            if (mm == N)//т.е. матрица сингулярная
                return false;
            double helpvar;
            for (uint c = 0; c < N; c++)
            {
                helpvar = A[m][c];
                A[m][c] = A[mm][c];
                A[mm][c] = helpvar;
            }
            helpvar = b[m];
            b[m] = b[mm];
            b[mm] = helpvar;
        }
        det *= A[m][m];
        // делим m строку на первый (диаг.) элемент
        // т.е. делаем диаг.элементы единицами
        double q = A[m][m];
        for (uint c = m; c < N; c++)
            A[m][c] /= q;
        b[m] /= q;
        // отнимаем от всех строк (кроме m-ой) m-ю строку умноженную на первый элемент текущей (r-ой) строки
        // т.е. обнуляем за исключением диаг. элемента m столбец
        for (uint r = 0; r < N; r++)//r номер строки
        {
            if (r == m)
                continue;
            q = A[r][m];
            for (uint c = m; c < N; c++)//c номер столбца
                A[r][c] -= q*A[m][c];//A'[r][c] = A[r][c] - A[r][m]*A[m][c];
            b[r] -= q*b[m];//b'[r] = b[r] - A[r][m]*b[m];
        }
        // отображение СЛАУ после текущей итерации
        // cout <<"m = "<<m<< endl;
        // showSLE(A,b,N);
    }
    det *= A[N-1][N-1];
    //cout <<"det = "<<det<< endl;
    if (det == 0.0 || std::isnan(det))//(fabs(det)<1e-5)
    {
        std::cout <<"The matrix is a singular!: det = "<<det<<std::endl;
        return false;
    }
    // находим решение
    for (uint r = 0; r < N; r++)
        x[r] = b[r];

    return true;
}
///
/// \details Решение системы линейных алгебраических уравнений (СЛАУ)
/// с использованием обраной матрицы.
/// \param[out] x - вектор корней - решение СЛАУ
/// \param[in] A - матрица
/// \param[in] b - вектор свободных членов
/// \param[in] N - размерность СЛАУ
/// \return bool - true/false - успешные/неуспешные вычисления
/// \warning Матрица A и вектор свободных членов b в процессе расчёта изменяются !
bool SLE_InverseMatrix(double* x, double** A, double* b, const uint N)
{
    double det;// детерминант
    double** iA = new double* [N];
    for (uint i = 0; i < N; i++)
        iA[i] = new double [N];

    bool res = getInverseMatrix(iA, A, N, det);
    if (res)
        getMatrixVectorProduct(x, iA, b, N);

    for (uint i = 0; i < N; i++)
        delete [] iA[i];
    delete [] iA;

    return res;
}
///
/// \details Вычисление обратной матрицы - основано на методе обращения
/// используемом в методе Гаусса - Жордана.
/// Решение СЛАУ с использованием обратной матрицы удобно, когда необходимо решать много СЛАУ
/// с одной и той же матрицей A, но с разными векторами свободных членов.
/// Основной недостаток этого метода (с использованием обратной матрицы) то, что он наиболее
/// медленный из всех методов.
/// \param[out] iA - обратная матрица
/// \param[in] A - матрица
/// \param[in] N - размерность СЛАУ
/// \return bool - true/false - успешные/неуспешные вычисления
/// \warning Матрица A и вектор свободных членов b в процессе расчёта изменяются !
///
bool getInverseMatrix(double** iA, double** A, const uint N, double& det)
{
    // подразумеваемый формат матрицы: A[row][column]
    // делаем матрицу iA предварительно единичной, после обращения
    // она трансформируется в обратную
    for (uint r = 0; r < N; r++)// r номер строки
        for (uint c = 0; c < N; c++)// c номер столбца
            if (r != c)
                iA[r][c] = 0.0;
            else
                iA[r][c] = 1.0;
    // обращаем матрицы: A*A^-1 = E; => E*A^-1 = A^-1;
    det = 1.0;
    for (uint m = 0; m < N; m++)// m номер опорной строки или диаг. элемента
    {
        // перестановка строк, если опорный элемент нулевой
        if (A[m][m] == 0.0)
        {
            uint mm = m;
            while (mm < N && A[mm][m] == 0.0)
                mm++;
            if (mm == N)// т.е. матрица сингулярная
            {
                det = 0.0;
                break;
            }
            double helpvar;
            for (uint c = 0; c < N; c++)
            {
                helpvar = A[m][c];
                A[m][c] = A[mm][c];
                A[mm][c] = helpvar;
                helpvar = iA[m][c];
                iA[m][c] = iA[mm][c];
                iA[mm][c] = helpvar;
            }
        }
        det *= A[m][m];
        // делим m строку на первый (диаг.) элемент
        // т.е. делаем диаг.элементы единицами
        double q = A[m][m];//вспомогательная переменная
        for (uint c = m; c < N; c++)
            A[m][c] /= q;
        for (uint c = 0; c < N; c++)
            iA[m][c] /= q;
        // отнимаем от всех строк (кроме m-ой) m-ю строку умноженную на первый элемент
        // текущей (r-ой) строки, т.е. обнуляем m-ый столбец за исключением диаг. элемента
        for (uint r = 0; r < N; r++)//r номер строки
        {
            if (r == m)
                continue;
            q = A[r][m];
            for (uint c = 0; c < N; c++)// c номер столбца
            {
                A[r][c] -= q*A[m][c];//A'[r][c] = A[r][c] - A[r][m]*A[m][c];
                iA[r][c] -= q*iA[m][c];
            }
        }
        // отображение трансформации ед. в обратную матрицу
        // на текущей итерации
        // cout <<"m = "<<m<< endl;
        // showMatrix(iA, N);
    }
    det *= A[N-1][N-1];
    if (det == 0.0 || std::isnan(det))
    {
        //cout <<"The matrix is a singular!: det = "<<det<< endl;
        return false;
    }

    return true;
}
///
/// матрично-векторное умножение: result = M*v
/// \param[out] result - результат умножения матрицы на вектор
/// \param[in] M - матрица
/// \param[in] v - вектор
/// \param[in] N - размерность матриц
///
void getMatrixVectorProduct(double* result, double** M, double* v, const uint N)
{
    // подразумеваемый формат матрицы: A[row][column]
    for (uint k = 0; k < N; k++)// номер строки
    {
        result[k] = 0.0;
        for (uint n = 0; n < N; n++)// номер столбца
            result[k] += M[k][n]*v[n];
    }
}
///
/// матричное умножение: C = A*B
/// \param[out] C - матрица - результат умножения первой на вторую
/// \param[in] A - первая матрица
/// \param[in] B - вторая матрица
/// \param[in] N - размерность матриц
///
void getMatrixProduct(double** C, double** A, double** B, const uint N)
{
    // подразумеваемый формат матрицы: A[row][column]
    for (uint k = 0; k < N; k++)// номер строки
        for (uint n = 0; n < N; n++)// номер столбца
        {
            C[k][n] = 0.0;
            for (uint i = 0; i < N; i++)// номер строки и столбца
                C[k][n] += A[k][i]*B[i][n];
        }
}
///
/// \details Вычислется нижняя треугольная матрица L из разложения Холецкого
/// симметричной (A[i,j] = A[j,i]) положительно-определенной матрицы  A = L*L'.
/// \param[out] L - нижняя треугольная матрица
/// \param[in] A - факторизуемая матрица
/// \param[in] N - размерность матриц
/// \return bool - true/false - успешные/неуспешные вычисления
///
bool Cholesky_Decomposition(double** L, double** A, const uint N)
{
    double det = 1.0;// детерминант L
    // подразумеваемый формат матрицы: A[row][column]
    for (uint i = 0; i < N; i++)// rows: 0..N-1
    {
        double S;
        // недиагональные элементы
        for (uint j = 0; j < i; j++)// columns: 0..i-1
        {
            S = 0.0;
            for (uint k = 0; k < j; k++)
                S += L[i][k]*L[j][k];
            L[i][j] = (A[i][j] - S) / L[j][j];
            L[j][i] = 0.0;
        }
        // диагональные элементы
        S = 0.0;
        for (uint k = 0; k < i; k++)
            S += L[i][k]*L[i][k];
        L[i][i] = sqrt(A[i][i] -  S);
        det *= L[i][i]*L[i][i];
        // тест на существование факторизации
        if (det == 0.0 || std::isnan(det))
        {
            //cout <<"матрица A не факторизуема !"<< endl;
            return false;
        }        
    }

    return true;
}
///
/// \details Решение системы линейных алгебраических уравнений (СЛАУ): \f[A\cdot x = b \f]
/// с симметричной положительно-определенной матрицей A с помощью факторизации Холецкого (A = L*L' -> L):
/// A*x=b \sim L*L'*x=b \sim L*y = b -> y; L'*x = y -> x.
/// Основное достоинство - высокая скорость вычислений, недостаток - применим только для вышеуказанных матриц.
/// \param[out] x - вектор корней - решение СЛАУ
/// \param[in] A - матрица
/// \param[in] b - вектор свободных членов
/// \param[in] N - размерность СЛАУ
/// \return bool - true/false - успешные/неуспешные вычисления
///
bool SLE_Cholesky(double* x, double** A, double* b, const uint N)
{
    bool result;
    // создаём нижнюю треугольную матрицу L
    double** L = new double* [N];
    for (uint i = 0; i < N; i++)
        L[i] = new double [N];
    // делаем разложение Холецкого
    bool res = Cholesky_Decomposition(L, A, N);
    if (res)// матрица L существует
    {
        // решаем промежуточную СЛАУ L*y = b; => y
        double* y = new double [N];
        for (uint i = 0; i < N; i++)
        {
            double S = 0;
            for (uint j = 0; j < i; j++)
                S += L[i][j]*y[j];
            y[i] = (b[i] - S) / L[i][i];
        }
        // решаем приведенную СЛАУ L'*x = y; => x
        for (int i = static_cast<int>(N)-1; i > -1; i--)
        {
            double S = 0.0;
            for (int j = static_cast<int>(N)-1; j > i; j--)
                S += L[j][i]*x[j];
            x[i] = (y[i] - S) / L[i][i];
        }
        delete [] y;
        result = true;
    }
    else// матрица L не существует
    {
        result = false;
    }
    // удаляем матрицу L
    for (uint i = 0; i < N; i++)
        delete[] L[i];
    delete[] L;

    return result;
}
///
/// \details Генератор симметричной заведомо положительно-определенной
/// матрицы для проверки алгоритма факторизации Холецкого
///
//void generateSymmPosDefiniteMatrix(double** L, double** A, const uint N)
void generateSymmPosDefiniteMatrix(double** A, const uint N)
{
    double**  L = new double* [N];
    for (uint i = 0; i < N; i++)
        L[i] = new double [N];

    // подразумеваемый формат матрицы: A[row][column]
    // рассчитываем нижнюю треугольную матрицу
    for (uint i = 0; i < N; i++)
        for (uint j = i; j < N; j++)
        {
            L[j][i] = ceil(Amp*rnd());
            if (i != j)
                L[i][j] = 0.0;// наддиагональная часть матрицы
            else
                if (L[i][j] == 0.0)// диагон. элемент - нулевой
                    L[i][j]++;// делаем ненулевым, иначе матрица не полож.-опред-ная
        }
    // рассчитываем A = L*L'
    for (uint i = 0; i < N; i++)
        for (uint j = 0; j < N; j++)
        {
            A[i][j] = 0.0;
            for (uint k = 0; k < N; k++)
                A[i][j] += L[i][k]*L[j][k];
        }
    // отображаем L для последующего сравнения с результатом расчета
    // функции Cholesky_Decomposition(A)
    //    cout <<"original lower triangular matrix L = "<< endl;
    //    showMatrix(L, N);
    for (uint i = 0; i < N; i++)
        delete [] L[i];
    delete [] L;
}
///
/// \details Отображение СЛАУ
///
void showSLE(double** A, double* b, const uint N)
{
    if (N > 8)
        return;
    for (uint i = 0; i < N+2; i++)
        std::cout <<SLINE;
    std::cout << std::endl;
    for (uint i = 0; i < N; i++)
    {
        for (uint j = 0; j < N; j++)
            std::cout <<A[i][j]<<"\t";
        std::cout <<"|\t"<<b[i]<< std::endl;
    }
    for (uint i = 0; i < N+2; i++)
        std::cout <<SLINE;
    std::cout << std::endl;
}
///
/// \details Отображение матрицы
///
void showMatrix(double** A, const uint N)
{
    if (N > 8)
        return;
    for (uint i = 0; i < N; i++)
        std::cout << SLINE;
    std::cout << std::endl;
    for (uint i = 0; i < N; i++)
    {
        for (uint j = 0; j < N; j++)
            std::cout <<A[i][j]<<"\t";
        std::cout << std::endl;
    }
    for (uint i = 0; i < N; i++)
        std::cout << SLINE;
    std::cout << std::endl;
}
///
/// \details Генератор случайной целочисленной матрицы
///
void generateRandMatrix(double** M, const uint N)
{
    for (uint i = 0; i < N; i++)
        for (uint j = 0; j < N; j++)
            M[i][j] = ceil(Amp*rnd());
}
///
/// \details Генератор случайной целочисленной симметричной матрицы
///
void generateSymmMatrix(double** M, const uint N)
{
    for (uint i = 0; i < N; i++)
        for (uint j = i; j < N; j++)
            M[i][j] = M[j][i] = ceil(Amp*rnd());
}
///
/// \details Генератор случайного целочисленного вектора
///
void generateRandVec(double* v, const uint N)
{
    for (uint i = 0; i < N; i++)
        v[i] = ceil(Amp*rnd());
}
///
// полиномиальное сглаживание минимизацией стандартного отклонения
/// \param [out] {xs, ys} - вектора сглаженных данных;
/// \param [in] Ns - размер векторов сглаженных данных;
/// \param [in] {x, y} - вектора  исходных данных;
/// \param [in] N - размер векторов исходных данных;
/// \param [in] n - порядок сглаживающего полинома:
///     - если n > 0, то n - максимально возможный порядок сглаживающего полинома:
///             \f$y(x)\ =\ a_n\cdot x^n + ... + a_1\cdot x + a_0\f$,
///       в процессе вычислений подбирается оптимальный порядок которого не больше n;
///     - если n < 0, то фиксированный порядок сглаживающего полинома равен |n|;
/// \return стандартное отклонение данных для найденного сглаживающего полинома
/// \warning
/// -# исходный вектор координат x следует сдвигать в начало координат,
///  иначе матричные элементы переполняются, что ведет к полной потере точности;
/// -# матрица коэф-тов для сглаживающих полиномов не факторизуема,
/// если порядок этих полиномов превышает 11 (от кол-ва данных этот факт не зависит);
///
double PolynomSmoothing(double* xs, double* ys, const uint Ns, double* x, double* y, const uint N, const int n)
{
    double STD = -1.0;// стандартное отклонение
    std::vector<double> a;// коэффициенты сглаживающего полинома
    if (n > 0)// случай оптимального сглаживающего полинома
    {
        // вычисляем полином. коэффициенты a и стандартное отклонение
        STD = getCoeffOptimalSmoothPoly(a, static_cast<uint>(n)+1, x, y, N);
        // cout <<"порядок оптим. сглаживающего полинома: n = "<<a.size()-1<<"\t STD = "<<STD<< endl;
        if (a.size() != 0)
        {
            double dxs = (x[N-1] - x[0]) / (Ns-1);
            for (uint j = 0; j < Ns; j++)
            {
                xs[j] = x[0] + j*dxs;
                ys[j] = Polynomial(xs[j], a);
//                cout <<j<<"\t"<<xs[j]<<"\t"<<ys[j]<< endl;
            }
        }
    }
    else// случай сглаживающего полинома фиксированного порядка
    {
        // задаем размерность вектора коэфф-тов сглаживающего полинома
        a.resize( static_cast<uint>(-n)+1 );
        // вычисляем полином. коэффициенты a и стандартное отклонение
        double dxs = (x[N-1] - x[0]) / (Ns-1);
        double dx = (x[N-1] - x[0]) / (N-1);
        int k = static_cast<int>(round(0.5*dxs/dx));// кол-во точек слева и справа
        if (k == 0)//нижний предел k
            k = 1;
        //cout <<"k = "<<k<< endl;
        double* xx = new double [2*k+1];
        double* yy = new double [2*k+1];
        int j, jl, jr;
        for (int i = 0; i < static_cast<int>(Ns); i++)
        {
            xs[i] = x[0] + i*dxs;
            // индекс ближайшей точки из вектора сглаживаемых данных к текущей
            // xs[i] = x[j]
            j = static_cast<int>(round(i*dxs/dx));
            jl = j - k; jr = j + k;
            if (jl < 0)
            {
                jl = 0;
                jr = 2*k;
            }
            if (jr > static_cast<int>(N)-1)
            {
                jl = static_cast<int>(N)-1 - 2*k;
                jr = static_cast<int>(N)-1;
            }
            // формируем массивы локальных данных {xx,yy}
            for (int jj = jl; jj < jr+1; jj++)
            {
                xx[jj-jl] = x[jj];
                yy[jj-jl] = y[jj];
            }
            // вычисляем полином. коэффициенты a и стандартное отклонение
            STD = getCoeffSmoothPoly(a, xx, yy, static_cast<uint>(2*k+1));
            ys[i] = Polynomial(xs[i], a);
        }
        delete[] xx;
        delete[] yy;
    }

    return STD;
}
///
/// \details Вычисляются коэффициенты сглаживающего полинома оптимального порядка,
/// но порядка не выше n, а также стандартное отклонение для N-мерного вектора данных
/// сглаженного с помощью этого полинома.
/// \param[out] a - вектор коэффициентов сглаживающего полинома
/// \param[in]  n - максимально возможный размер вектора коэффициентов a,
///             т.е. n-1 максимально возможный порядок сглаживающего полинома
/// \param[in]  {x,y} - вектора данных
/// \param[in]  N - размерность векторов данных
/// \return стандартное отклонение данных для вычисленного сглаживающего полинома
///
double getCoeffOptimalSmoothPoly(std::vector<double>& a, uint n, const double* x, const double* y, const uint N)
{
    // создаём и вычисляем матрицу коэффициентов X[k][j]
    // при коэффициентах a[j] сглаживающего полинома
    double** X = new double* [n];
    for (uint i = 0; i < n; i++)
        X[i] = new double [n];
    // и вектор свободных членов f[k] = Sum(y[i]*x[i]^k)
    double* f = new double [n];
    getCoeffMatrix(X, f, n, x, y, N);
    // cout <<"X ="<< endl;
    // showSLE(X, f, n);
    // производим поиск оптимального сглаживающего полинома
    // находим коэффициенты, решая СЛАУ с помощью
    // разложения Холецкого, т.к. получившаяся матрица
    // симметричная (всегда ли положительно определенная ???),
    // опираясь на предварительно вычисленные матрицу X и вектор f
    std::vector<double> ca(1);// коэфф. сглаживающего полинома текущего l-го порядка
    double mean_y = f[0]/N;// среднее вектора данных у
    double STD = 0.0;// среднеквадратичное отклонение
    for (uint i = 0; i < N; i++)
        STD += (y[i] - mean_y)*(y[i] - mean_y);
    ca[0] = 0;  // нулевой порядок сглаживающего полинома
    double** hX;// вспомогательная матрица коэфф-тов для текущего порядка сглаж.полинома
    double* hf; // вспомогательный вектор своб.членов
    double* coeff;
    double cSTD;// текущее среднеквадратичное отклонение
    for (uint l = 1; l < n; l++)// l - текущий порядок полинома
    {
        hX = new double* [l];
        for (uint i = 0; i < l; i++)
            hX[i] = new double [l];
        hf = new double [l];
        // заполняем матрицу коэффициентов и вектор своб.членов
        // текущей размерности
        for (uint i = 0; i < l; i++)
        {
            hf[i] = f[i];
            for (uint j = i; j < l; j++)
                hX[i][j] = hX[j][i] = X[i][j];
        }
        // находим решение СЛАУ - коэффициенты сглаживающего
        // полинома текущего l-го порядка
        coeff = new double [l];
        bool res = SLE_Cholesky(coeff, hX, hf, l);
        if (res)// коэффициенты успешно вычислены
        {
            if (ca.size() != 0)
                ca.resize(l);
            for (uint i = 0; i < l; i++)
                ca[i] = coeff[i];
            cSTD = 0;
            double Y;
            for (uint i = 0; i < N; i++)
            {
                Y = Polynomial(x[i], ca);
                cSTD += (Y - y[i])*(Y - y[i]);
            }
            if (STD > cSTD)
            {
                STD = cSTD;
                a = ca;// новый вектор коэф-тов оптимального сглаж. полинома
            }
            //cout <<l-1<<" : "<<cSTD<< endl;
        }
        else// прерываем цикл - если текущая не факторизуема, ??? так ли
        {// то и последующие матрицы тоже не факторизуемы
            delete[] coeff;
            break;
        }
        delete[] coeff;
        for (uint i = 0; i < l; i++)
            delete[] hX[i];
        delete[] hX;
        delete[] hf;
    }
    // уничтожаем массивы X и f
    for (uint i = 0; i < n; i++)
        delete[] X[i];
    delete[] X;
    delete[] f;

    return STD;
}
///
/// \details Вычисляются коэффициенты сглаживающего полинома порядка n = a.size()-1
/// Y(x) = a[0]+a[1]*x+a[2]*x^2+...a[n]*x^n, а также стандартное отклонение
/// для N-мерного вектора данных сглаженного с помощью этого полинома.
/// \param[out] a - вектор коэффициентов сглаживающего полинома
///                 размер a должен быть задан ! - он определяет порядок аппроксим. полинома
///                 если нужна сшивка на левом конце интервала с предыдущим сглаживающим полиномом,
///                 то д.б. заданы a[0] и a[1] текущего сглаживающего полинома, используя предыдущий полином в точке x = x0
/// \param[in]  {x,y} - вектора данных
/// \param[in]  N - размерность векторов данных
/// \return стандартное отклонение данных для вычисленного сглаживающего полинома
/// https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%BD%D0%B0%D0%B8%D0%BC%D0%B5%D0%BD%D1%8C%D1%88%D0%B8%D1%85_%D0%BA%D0%B2%D0%B0%D0%B4%D1%80%D0%B0%D1%82%D0%BE%D0%B2
///
double getCoeffSmoothPoly(std::vector<double>& a, const double* x, const double* y, const uint N)
{
    // создаём и вычисляем матрицу коэффициентов X[k][j]
    // при коэффициентах a[j] сглаживающего полинома
    uint n = static_cast<uint>(a.size());// число коэфф-тов сглаживающего полинома
    double** X = new double* [n];
    for (uint i = 0; i < n; i++)
        X[i] = new double [n];
    // и вектор свободных членов f[k] = Sum(y[i]*x[i]^k)
    double* f = new double [n];
    getCoeffMatrix(X, f, n, x, y, N);
    /// особый случай: при необходимости удовлетворения условиям сшивки текущего полинома
    /// на одной из границ:
    /// y(x0) = Poly(x = x0, a) = Sum_{n=0}^m { a[n]*(x - x0)^n } = a[0]
    /// y'(x0) = Poly'(x = x0, a) = Sum_{n=1}^m { a[n]*n*(x - x0)^(n-1) } = a[1]
    /// трансформируем СЛАУ так, чтобы воспроизводились заданные a[0] и a[1]
    if (a.size() > 2 && a.at(0) != 0 && a.at(1) != 0)
    {
        /// первые 2 строки и столбца становятся диагональными и единичными
        for (uint i = 0; i < 2; i++)
            for (uint j = 0; j < n; j++)
                X[i][j] = X[j][i] = (i == j) ? 1.0 : 0.0;
        /// теперь коррекция свободных членов
        f[0] = a[0];
        f[1] = a[1];
        double* xm = new double [N];
        for (uint l = 0; l < 2; l++)
            for (uint k = 0; k < N; k++)
                xm[k] = (l > 0) ? xm[k]*x[k] : 1.0;
        for (uint l = 2; l < n; l++)
            for (uint k = 0; k < N; k++)
            {
                xm[k] = (l > 0) ? xm[k]*x[k] : 1.0;
                f[l] -= (a[0] + a[1]*x[k])*xm[k];
            }
        delete [] xm;
    }
    // cout <<"X ="<< endl;
    // showSLE(X, f, n);
    // решаем СЛАУ методом Холецкого и сохраняем найденные коэфф-ты
    // сглаживающего полинома в вектор a
    double* coeff = new double [n];
    SLE_Cholesky(coeff, X, f, n);
    for (uint i = 0; i < n; i++)
        a[i] = coeff[i];
    // освобождаем память
    delete[] coeff;
    for (uint i = 0; i < n; i++)
        delete[] X[i];
    delete[] X;
    delete[] f;
    // находим стандартное отклонение
    double STD = 0.0, Y;
    for (uint i = 0; i < N; i++)
    {
        Y = Polynomial(x[i], a);
        STD += (Y - y[i])*(Y - y[i]);
    }

    return STD;
}
double getCoeffSmoothPoly(std::vector<double>& a, const std::vector<QPointF> points)
{
    const uint N = points.size();
    // создаём и вычисляем матрицу коэффициентов X[k][j]
    // при коэффициентах a[j] сглаживающего полинома
    uint n = static_cast<uint>(a.size());// число коэфф-тов сглаживающего полинома
    double** X = new double* [n];
    for (uint i = 0; i < n; i++)
        X[i] = new double [n];
    // и вектор свободных членов f[k] = Sum(y[i]*x[i]^k)
    double* f = new double [n];
    getCoeffMatrix(X, f, n, points);
    /// особый случай: при необходимости удовлетворения условиям сшивки текущего полинома
    /// на одной из границ:
    /// y(x0) = Poly(x = x0, a) = Sum_{n=0}^m { a[n]*(x - x0)^n } = a[0]
    /// y'(x0) = Poly'(x = x0, a) = Sum_{n=1}^m { a[n]*n*(x - x0)^(n-1) } = a[1]
    /// трансформируем СЛАУ так, чтобы воспроизводились заданные a[0] и a[1]
    if (a.size() > 2 && a.at(0) != 0 && a.at(1) != 0)
    {
        /// первые 2 строки и столбца становятся диагональными и единичными
        for (uint i = 0; i < 2; i++)
            for (uint j = 0; j < n; j++)
                X[i][j] = X[j][i] = (i == j) ? 1.0 : 0.0;
        /// теперь коррекция свободных членов
        f[0] = a[0];
        f[1] = a[1];
        double* xm = new double [N];
        for (uint l = 0; l < 2; l++)
            for (uint k = 0; k < N; k++)
                xm[k] = (l > 0) ? xm[k]*points[k].x() : 1.0;
        for (uint l = 2; l < n; l++)
            for (uint k = 0; k < N; k++)
            {
                xm[k] = (l > 0) ? xm[k]*points[k].x() : 1.0;
                f[l] -= (a[0] + a[1]*points[k].x())*xm[k];
            }
        delete [] xm;
    }
    // cout <<"X ="<< endl;
    // showSLE(X, f, n);
    // решаем СЛАУ методом Холецкого и сохраняем найденные коэфф-ты
    // сглаживающего полинома в вектор a
    double* coeff = new double [n];
    SLE_Cholesky(coeff, X, f, n);
    for (uint i = 0; i < n; i++)
        a[i] = coeff[i];
    // освобождаем память
    delete[] coeff;
    for (uint i = 0; i < n; i++)
        delete[] X[i];
    delete[] X;
    delete[] f;
    // находим стандартное отклонение
    double STD = 0.0, Y;
    for (uint i = 0; i < N; i++)
    {
        Y = Polynomial(points[i].x(), a);
        STD += (Y - points[i].y())*(Y - points[i].y());
    }

    return STD;
}
/// расчёт матрицы X и вектора f, используемых для нахождения коэффициентов
/// сглаживающего полинома при решении СЛАУ:
///         X*a = f { X[k][j]*a[j] = f[k] = Sum(y[i]*x[i]^k) }
/// n - размерность матрицы X и длина вектора f
///
void getCoeffMatrix(double** X, double* f, const uint n, const double* x, const double* y, const uint N)
{
    // вспомогательный вектор, содержащий вектор
    // координат x в m(=j+k)-ой степени: xm = x[i]^m
    double* xm = new double [N];
    double Sxm;
    // расчёт элементов нулевой строки матрицы X: X[k=0][j]
    // начальная предустановка x[i]^0 = 1
    f[0] = 0.0;
    for (uint i = 0; i < N; i++)
    {
        xm[i] = 1;
        f[0] += y[i];
    }
    X[0][0] = N;
    for (uint j = 1; j < n; j++)
    {
        f[j] = 0.0;
        for (uint i = 0; i < N; i++)
        {
            xm[i] *= x[i];
            f[j] += y[i]*xm[i];
        }
        Sxm = 0;
        for (uint i = 0; i < N; i++)
            Sxm += xm[i];
        X[0][j] = Sxm;
    }
    // расчёт элементов k-ой строки матрицы X: X[k>0][j]
    for (uint k = 1; k < n; k++)
    {
        // каждая строка совпадает с предыдущей,
        // сдвинутой влево на один элемент
        for (uint j = 0; j < n-1; j++)//j = 0..n-2
            X[k][j] = X[k-1][j+1];
        // отдельно рассчитаем последний элемент строки X[k][n]
        for (uint i = 0; i < N; i++)
            xm[i] *= x[i];
        Sxm = 0;
        for (uint i = 0; i < N; i++)
            Sxm += xm[i];
        X[k][n-1] = Sxm;
    }
    delete[] xm;
}
void getCoeffMatrix(double** X, double* f, const uint n, const std::vector<QPointF>& data)
{
    // вспомогательный вектор, содержащий вектор
    // координат x в m(=j+k)-ой степени: xm = x[i]^m
    double* xm = new double [data.size()];
    double Sxm;
    // расчёт элементов нулевой строки матрицы X: X[k=0][j]
    // начальная предустановка x[i]^0 = 1
    f[0] = 0.0;
    for (size_t i = 0; i < data.size(); i++)
    {
        xm[i] = 1;
        f[0] += data[i].y();
    }
    X[0][0] = data.size();
    for (uint j = 1; j < n; j++)
    {
        f[j] = 0.0;
        for (size_t i = 0; i < data.size(); i++)
        {
            xm[i] *= data[i].x();
            f[j] += data[i].y()*xm[i];
        }
        Sxm = 0;
        for (size_t i = 0; i < data.size(); i++)
            Sxm += xm[i];
        X[0][j] = Sxm;
    }
    // расчёт элементов k-ой строки матрицы X: X[k>0][j]
    for (uint k = 1; k < n; k++)
    {
        // каждая строка совпадает с предыдущей,
        // сдвинутой влево на один элемент
        for (uint j = 0; j < n-1; j++)//j = 0..n-2
            X[k][j] = X[k-1][j+1];
        // отдельно рассчитаем последний элемент строки X[k][n]
        for (size_t i = 0; i < data.size(); i++)
            xm[i] *= data[i].x();
        Sxm = 0;
        for (size_t i = 0; i < data.size(); i++)
            Sxm += xm[i];
        X[k][n-1] = Sxm;
    }
    delete[] xm;
}
// возвращает значение полинома (a.size()-1)-го порядка с вектором коэффициентов a
double Polynomial(double x, const std::vector<double>& a)
{
//    double y = 0.0, xn = 1.0;
//    for (uint n = 0; n < a.size(); n++)
//    {
//        y += a[n]*xn;
//        xn *= x;
//    }
    double y = a[a.size()-1];
    for (int n = static_cast<int>(a.size())-2; n > -1; n--)
        y = y*x + a[static_cast<uint>(n)];

    return y;
}
// возвращает значение полинома n-го порядка с вектором корней roots
double PolynomialR(double x, const std::vector<double>& roots)
{
    double y = 1.0;
    for (uint n = 0; n < roots.size(); n++)
        y *= (x - roots[n]);
    return y;
}
///
/// \details Запись в файл
///
void WriteFile(const char* FileName, double* x, double* y, const int N)
{
    std::ofstream oFile(FileName);

    std::cout <<"\nзапись в файл: "<<FileName<< std::endl;
    for (int i = 0; i < N; i++)
    {
        oFile<<x[i]<<"\t"<<y[i]<< std::endl;
        //cout <<i<<" : "<<x[i]<<"\t"<<y[i]<< endl;
    }
    oFile.close();
}
// чтение из файла ???
// фактически функция - не рабочая т.к. массивы x,y
// после возврата из ф-ции - становятся не действительными
void ReadFile(const char* FileName, double* x, double* y, int& N)
{
    std::ifstream iFile(FileName);
    double xx,yy;
    std::ios::pos_type pos;

    std::cout <<"\nчтение из файла: "<<FileName<< std::endl;
    N = 0;
    while(!iFile.eof())
    {
        N++;
        iFile>>xx>>yy;
        pos = iFile.tellg();
        std::cout <<N<<" : "<<pos<<" : "<<xx<<"\t"<<yy<<"\t"<< std::endl;
    }
    N--;
    std::cout <<"N = "<<N<< std::endl;
    //передвигаем указатель на начальную позицию
    //iFile.seekg(0, ios_base::beg);//не работает
    //iFile.seekg(start_pos);//не работает
    iFile.close();
    iFile.open(FileName, std::ios_base::in);
    x = new double [N];
    y = new double [N];
    for (int i = 0; i < N; i++)
    {
        iFile>>x[i]>>y[i];
        std::cout <<i<<" : "<<x[i]<<"\t"<<y[i]<< std::endl;
    }
    iFile.close();
}
///
/// \details Запись в файл
///
void WriteFile(const char* FileName, const std::vector<double>& x, const std::vector<double>& y)
{
    std::ofstream oFile(FileName);

    std::cout <<"\nзапись в файл: "<<FileName<< std::endl;
    for (uint i = 0; i < y.size(); i++)
        oFile<<x[i]<<"\t"<<y[i]<< std::endl;
    oFile.close();
}
///
/// \details Чтение из файла
///
void ReadFile(const char* FileName, std::vector<double>& x, std::vector<double>& y)
{
    std::ifstream iFile(FileName);
    double xx,yy;
    ulong N = 0;
    // определяем кол-во строк в файле, т.е. длину векторов
    std::cout <<"\nчтение из файла: "<<FileName<< std::endl;
    while (!iFile.eof())
    {
        N++;
        iFile>>xx>>yy;
    }
    N--;
    // переопределяем длину векторов
    x.resize(N);
    y.resize(N);
    // передвигаем указатель на начальную позицию
    iFile.close();
    iFile.open(FileName, std::ios_base::in);// правильные варианты не работают - см.выше
    // заполняем вектора данными из файла
    for (uint i = 0; i < N; i++)
        iFile>>x[i]>>y[i];
    iFile.close();
}
///
/// \details Запись в файл данных из стандартного вектора типа QPointF
///
//void WriteFile(const char* FileName, const std::vector<QPointF>& data)
//{
//    std::ofstream oFile(FileName);

//    std::cout <<"\nзапись в файл: "<<FileName<< std::endl;
//    for (uint i = 0; i < data.size(); i++)
//        oFile<<data[i].x()<<"\t"<<data[i].y()<< std::endl;
//    oFile.close();
//}
///
/// \details Чтение из файла в стандартный вектор данных типа QPointF
///
//void ReadFile(const char* FileName, std::vector<QPointF>& data)
//{
//    std::ifstream iFile(FileName);
//    double xx,yy;
//    ulong N = 0;
//    //определяем кол-во строк в файле - т.е. длину векторов
//    std::cout <<"\nчтение из файла: "<<FileName<< std::endl;
//    while (!iFile.eof())
//    {
//        N++;
//        iFile>>xx>>yy;
//    }
//    N--;
//    //переопределяем длину векторов
//    data.resize(N);
//    //передвигаем указатель на начальную позицию
//    iFile.close(); iFile.open(FileName, std::ios_base::in);//правильные варианты не работают - см.выше
//    //заполняем вектора данными из файла
//    for (uint i = 0; i < N; i++)
//    {
//        iFile>>xx>>yy;
//        data[i] = QPointF(xx, yy);
//    }
//    iFile.close();
//}
///
/// \details Факторизация числа числа f, делители которого сохраняются в вектор Factor
/// \param[in] f - факторизуемое число
/// \param[out] Factor - вектор множителей
///
void Factorization(long long f, std::vector<int>& Factor)
{
    long long int remainder = f;// текущий остаток от числа f
    long long int divider = 2;// текущий делитель числа f
    long long int rem = 1;// текущий сокращенный множитель (который уже отделен от f)
    ///
    if (remainder < 0)
    {
        Factor.push_back( -1 );
        remainder *=-1;
    }
    while (remainder != 1)
    {
        if (remainder % divider == 0) // проверяем - делится ли текущий остаток на текущий делитель
        {
            remainder /= divider;// сокращаем число на текущий делитель
            rem *= divider;// добавляем текущий делитель в сокращенный множитель
            Factor.push_back( static_cast<int>(divider) );
        }
        else// текущий остаток уже не делится на текущий делитель,
            // поэтому увеличиваем делитель до следующего числа, не являющегося произведением уже известных делителей
        {
            divider++;
            // пропускаем делители, присутствующие уже в сокращенном множителе
            while (rem % divider == 0)
                divider++;
        }
    }
}
///
/// \details Вычисление НОД на основе массивов их делителей
/// \param[in] Factor1 - вектор множителей первого числа
/// \param[in] Factor2 - вектор множителей второго числа
/// \return НОД двух чисел
/// \warning Содержимое вектора Factor2 меняется в процессе вычислений !
///
long long int getNOD(std::vector<int>& Factor1, std::vector<int>& Factor2)
{
    long long int NOD = 1;
    ///
    for (unsigned int i = 0; i < Factor1.size(); i++)
        for (unsigned int j = 0; j < Factor2.size(); j++)
            if (Factor1[i] == Factor2[j])
            {
                NOD *= Factor2[j];
                Factor2[j] = -1;//исключааем возможность повторного учёта этого множителя
                break;
            }
    return NOD;
}
///
/// \details Вычисление НОД двух чисел
/// \param[in] f1 - первое число
/// \param[in] f2 - второе число
/// \return НОД двух чисел
/// \warning Содержимое вектора Factor2 меняется в процессе вычислений !
///
long long int getNOD(long long int f1, long long int f2)
{
    std::vector<int> Factor1, Factor2;
    Factorization(f1, Factor1);
    Factorization(f2, Factor2);
    return getNOD(Factor1, Factor2);
}
