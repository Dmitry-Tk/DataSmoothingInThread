#include <vector>
#include <iostream>
#include <cmath>
#include "cfilter.h"
#include "mathlib.h"

using namespace std;


unsigned int CFilter::Nfilters = 0;/// число фильтров

///
// конструктор
/// \param[in] m       - порядок фильтра
/// \param[in] length  - длина фильтра - кол-во точек используемое для сглаживания
/// \param[in] N       - число каскадов фильтра
/// \param[in] FilterName - имя фильтра
///
CFilter::CFilter(unsigned int m, unsigned int length, unsigned int N, QString FilterName)
{
    genNewFilter( m, length, N, FilterName);
    Nfilters++;
    /// ключ: отображать ли состояние объекта в процессе работы (для отладки)
    ShowObjectInfo = false;
}
///
/// генератор нового фильтра
void CFilter::genNewFilter(unsigned int m, unsigned int length, unsigned int N, QString FilterName)
{
    /// имя фильтра
    if (FilterName == "")
        this->FilterName = "Poly: "+QString::number(m)+
                        "-"+QString::number(length)+"-"+QString::number(N);
//        this->FilterName = QString("filter #") + QString::number(Nfilters);
    else
        this->FilterName = FilterName;
    /// кол-во коэфф-тов в аппрокс. полиноме = порядок аппроксим. полинома-1
    this->m = m + 1;
    /// длина фильтра - число точек по которым вычисляется отфильтров. значение
    this->length = static_cast<int>(length);
    /// глубина фильтрации - число уровней
    this->N = N;
    /// производные параметры
    halflength = static_cast<unsigned int>(length)/2;// половина длины фильтра
    shift = N*halflength;// (- 1 перепроверить!!!) сдвиг между входным и фильтруемым отсчётом
    index = cindex = -1;// сброс счётчикв прокачанных данных и текущего индекса
    if (N == 1)
    {
        if (In != nullptr)
        {
            delete[] In;
            In = nullptr;
        }
        if (in != nullptr)
            delete[] in;
        in = new double[length];
    }
    else
    {
        if (in != nullptr)
            delete[] in;
        in = new double[N*static_cast<unsigned int>(length)];
        if (In != nullptr)
        {
            for (unsigned int i = 0; i < N; i++)
                delete[] In[i];
            delete[] In;
        }
        In = new double*[N];
        for (unsigned int i = 0; i < N; i++)
            In[i] = new double[length];
    }
    if (coeff != nullptr)
        delete[] coeff;
    coeff = new double[length];
    getFLoatCoeff();
    thSTD = 0.0;
    for (int i = 0; i < length; i++)
        thSTD += coeff[i]*coeff[i];
}
///
///
CFilter::~CFilter()
{
    if (in != nullptr)
    {
        delete[] in;
        in = nullptr;
    }
    if (N > 1)
    {
        if (In != nullptr)
        {
            for (unsigned int i = 0; i < N; i++)
                delete[] In[i];
            delete[] In;
            In = nullptr;
        }
    }
    if (coeff != nullptr)
    {
        delete[] coeff;
        coeff = nullptr;
    }
    if (ShowObjectInfo)
        cout << "The CFilter has been destructed !" << endl;
    Nfilters--;
}
///
/// ввод нового отсчёта в очередь
void CFilter::push(double newSample)
{
    index++;
    if (N == 1)
    {
        cindex++;
        if (cindex == length)
            cindex = 0;
        in[cindex] = newSample;
    }
    else
    {
        in[index % (static_cast<int>(N)*length)] = newSample;
        cindex++;
        if (cindex == length)
            cindex = 0;
        for (unsigned int n = 0; n < N; n++)
        {
            if (n == 0)
                In[n][cindex] = newSample;
            else
                In[n][cindex] = _getConvolution(n-1);// свёртка предыдущего уровня
        }
    }
}
///
// возвращение текущего последнего отсчёта из очереди фильтра
/// \return текущий последний отсчёт из очереди фильтра
///
double CFilter::GetLastSampleFromQueue()
{
    return in[shift];
}
///
// вычисление свёртки одно(много)каскадного фильтра
/// \return значение свёртки фильтра (текущий отфильтрованный отсчёт)
/// - для однокаскадного фильтра \f$ \sum\limits_{i=0}^{length-1}coeff_i\cdot in_{length-1-i}\f$
/// - для многокаскадного фильтра \f$ \sum\limits_{i=0}^{length-1}coeff_i\cdot In_{N-1,length-1-i}\f$
/// \warning Учитывая, что заполнение фильтра данными происходит конвейерным образом, то вышеприведенные
/// формулы трансформируются:
/// - для однокаскадного фильтра \f[\sum\limits_{j=0}^{cindex}coeff_{length-1-cindex-j}\cdot in_j +
///                                 \sum\limits_{j=cindex+1}^{length-1}coeff_{j-cindex-1}\cdot in_j\f]
/// - для многокаскадного фильтра \f[\sum\limits_{j=0}^{cindex}coeff_{length-1-cindex-j}\cdot In_{N-1,j} +
///                                  \sum\limits_{j=cindex+1}^{length-1}coeff_{j-cindex-1}\cdot In_{N-1,j}\f]
///
double CFilter::getConvolution()
{
    double conv = 0;

    if (N == 1)
    {
        for (int j = 0; j < cindex+1; j++)
            conv += coeff[length-1-cindex+j]*in[j];
        for (int j = cindex+1; j < length; j++)
            conv += coeff[j-cindex-1]*in[j];
    }
    else
    {
        for (int j = 0; j < cindex+1; j++)
            conv += coeff[length-1-cindex+j]*In[N-1][j];
        for (int j = cindex+1; j < length; j++)
            conv += coeff[j-cindex-1]*In[N-1][j];
    }

    return conv;
}
///
// вычисление свёртки на n-ом уровне
/// \return значение свёртки \f$n\f$-го каскада
/// \f$ \sum\limits_{i=0}^{length-1}coeff_i\cdot In_{n,length-1-i}\f$
/// \warning Учитывая, что заполнение фильтра данными происходит конвейерным образом, то вышеприведенная
/// формула трансформируется в \f[\sum\limits_{j=0}^{cindex}coeff_{length-1-cindex-j}\cdot In_{n,j} +
///                               \sum\limits_{j=cindex+1}^{length-1}coeff_{j-cindex-1}\cdot In_{n,j}\f]
///
double CFilter::_getConvolution(unsigned int n)
{
    double conv = 0;

    for (int j = 0; j < cindex+1; j++)
        conv += coeff[length-1-cindex+j]*In[n][j];
    for (int j = cindex+1; j < length; j++)
        conv += coeff[j-cindex-1]*In[n][j];

    return conv;
}
//
// очистка очереди
//
void CFilter::reset()
{
    for (int i = 0; i < length; i++)
        if (N == 1)
            in[i] = 0;
        else
            for (unsigned int il = 0; il < N; il++)
                In[il][i] = 0;
}
/* прежняя версия для целочисленного представления фильтра
 * в виде набора целых коэффициентов и целого знаменателя;
 * работает только для полиномов малых порядков и малых
 * значений числа точек
 * работающие варианты (вне них - фильтр не рассчитывается)
 * {_m+1: length(=2*n+1)} = {2: 5..11}; {3: 5..11,15..23}; {4: 7..11}; {5: 7}
void CFilter::GetCoeff()
{
    // создаём и вычисляем главную матрицу
    double** A = new double*[m];
    for (unsigned int i = 0; i < m; i++)
        A[i] = new double[m];
    std::vector<double> vec(2*(m-1));
    int _halflength = static_cast<int>(halflength);
//    qDebug()<<_halflength;
    for (unsigned int n = 0; n <= 2*(m-1); n++)
    {
        vec[n] = 0;
        for (int k = -_halflength; k <= _halflength; k++)
            vec[n] += pow(k, n);
//        qDebug()<<n<<"  "<<vec[n];
    }
    for (unsigned int r = 0; r < m; r++)
        for (unsigned int c = r; c < m; c++)
        {
            A[r][c] = vec[2*m-2 - r - c];
            A[c][r] = A[r][c];
        }
    // создаём и вычисляем вектор миноров
    double* M = new double[m];
    // создаём матрицу для алгебраических дополнений
    double** MA = new double*[m-1];
    for (unsigned int i = 0; i < m-1; i++)
        MA[i] = new double[m-1];
    for (unsigned int i = 0; i < m; i++)
    {   // заполняем вспомог. матрицу текущего минора
        for (unsigned int r = 0; r < m; r++)
            for (unsigned int c = 0; c < m-1; c++)
                if (r < i)
                    MA[r][c] = A[r][c];
                else
                    if (r == i)
                        continue;
                    else
                        MA[r-1][c] = A[r][c];
        M[i] = pow(-1.0, i+1+m)*GetDeterminant(MA, m-1);
    }
    // инициализируем вектор коэффициентов фильтра
    long long int* _coeff = new long long int[length];
    for (unsigned int i = 0; i < length; i++)
        _coeff[i] = 0;
    // вычисляем знаменатель фильтра
    long long int _denom = static_cast<long long int>(round(GetDeterminant(A, m)));
    // вычисляем вектор коэффициентов фильтра
    for (unsigned int n = 0; n < m; n++)
        for (int k = -_halflength; k <= _halflength; k++)
            _coeff[k+_halflength] += static_cast<long long int>(round(M[n]*pow(k,m-1-n)));
    // упрощение - сокращение множителей коэффициентов фильтра
    long long int NOD;
    long long int* _den = new long long int[length];
    //
    for (unsigned int i = 0; i < length; i++)
    {
        NOD = GetNOD(_coeff[i], _denom);// вычисляем НОД числителя и знаменателя текущего коэфф.
        _coeff[i] = _coeff[i] / NOD;// сокращаем на НОД числитель текущего коэфф.
        _den[i] = _denom / NOD;// сокращаем на НОД знаменатель текущего коэфф.
    }
    // вычисляем НОК всех знаменателей дробных коэфф. фильтра (знаменатель фильтра)
    // последовательно перебирая эти знаменатели
    long long int NOK = _den[0];
    for (unsigned int i = 1; i < length; i++)
        NOK = _den[i]*NOK / GetNOD(NOK, _den[i]);
    denom = static_cast<int>(NOK);
    // приводим к общему знаменателю все коэфф. фильтра
    for (unsigned int i = 0; i < length; i++)
        coeff[i] = denom/static_cast<int>(_den[i])*static_cast<int>(_coeff[i]);
    // очищаем память
    delete[] _den;
    delete[] _coeff;
    for (unsigned int i = 0; i < m-1; i++)
        delete [] MA[i];
    delete[] MA;
    for (unsigned int i = 0; i < m; i++)
        delete []A[i];
    delete[] A;
    delete[] M;
}
*/
///
/// \details здесь д.б. описание алгоритма генерации коэффициентов фильтра
///
void CFilter::getFLoatCoeff()
{
    /// создаём и вычисляем главную матрицу
    double** A = new double*[m];
    for (unsigned int i = 0; i < m; i++)
        A[i] = new double[m];
    std::vector<double> vec(2*(m-1));
    int _halflength = static_cast<int>(halflength);
    for (unsigned int n = 0; n <= 2*(m-1); n++)
    {
        vec[n] = 0;
        for (int k = -_halflength; k <= _halflength; k++)
            vec[n] += pow(k, n);
    }
    for (unsigned int r = 0; r < m; r++)
        for (unsigned int c = r; c < m; c++)
        {
            A[r][c] = vec[2*m-2 - r - c];
            A[c][r] = A[r][c];
        }
    /// создаём и вычисляем вектор миноров
    double* M = new double[m];
    /// создаём матрицу для алгебраических дополнений
    double** MA = new double*[m-1];
    for (unsigned int i = 0; i < m-1; i++)
        MA[i] = new double[m-1];
    for (unsigned int i = 0; i < m; i++)
    {   /// заполняем вспомог. матрицу текущего минора
        for (unsigned int r = 0; r < m; r++)
            for (unsigned int c = 0; c < m-1; c++)
                if (r < i)
                    MA[r][c] = A[r][c];
                else
                    if (r == i)
                        continue;
                    else
                        MA[r-1][c] = A[r][c];
        M[i] = pow(-1.0, i+1+m)*getDeterminant(MA, m-1);
    }
    /// инициализируем вектор коэффициентов фильтра
    double* _coeff = new double[length];
    for (int i = 0; i < length; i++)
        _coeff[i] = 0;
    /// вычисляем знаменатель фильтра
    double _denom = getDeterminant(A, m);
    /// вычисляем вектор коэффициентов фильтра
    for (unsigned int n = 0; n < m; n++)
        for (int k = -_halflength; k <= _halflength; k++)
            _coeff[k+_halflength] += M[n]*pow(k, m-1-n);
    /// приводим к общему знаменателю все коэфф. фильтра
    for (int i = 0; i < length; i++)
        coeff[i] = _coeff[i] / _denom;
    /// очищаем память
    delete[] _coeff;
    for (unsigned int i = 0; i < m-1; i++)
        delete[] MA[i];
    delete[] MA;
    for (unsigned int i = 0; i < m; i++)
        delete[] A[i];
    delete[] A;
    delete[] M;
}
///
// информация о фильтре
///
void CFilter::showFilterInfo(bool in_detail)
{
    if (!in_detail)
    {
        cout  << FilterName.toStdString();
        return;
    }
    cout  << "  Полиномиальный фильтр (CFilter):"
          << "\n    FilterName: " << FilterName.toStdString()
          << "\n    m = " << m-1 << " - порядок аппроксим. полинома"
          << "\n    length = " << length << " - число точек по которым происходит сглаживание"
          << "\n    levels = " << N << " - число каскадов в фильтре"
          << "\n    semilength = " << halflength
          << "\n    shift = " << shift << " - сдвиг (вносимый фильтром) между текущим и фильтров. отсчётами"
          << "\n    coeff = {" << endl;
    int _halflength = static_cast<int>(halflength);
    for (int k = -_halflength; k <=_halflength; k++)
        cout << "\t" << k << "\t" << coeff[k+_halflength] << endl;
    cout << "  \t}" << endl;
    cout << "\n    thSTD = " << thSTD << endl;
}
///
/// получить параметры фильтра
void CFilter::getFilterParams(unsigned int& m, unsigned int& length, unsigned int& N)
{
    /// возвращаем порядок полинома (!), а не кол-во коэфф-тов, также как и в конструкторе
    m = this->m-1;
    length = this->length;
    N = this->N;
}
///
// выполнение фильтрации вектора данных - вычисление свёртки фильтра и входного вектора данных
// краевые эффекты рассогласования произодных нивелированы
/// \param[in] inData   - входной вектор данных;
/// \param[out] outData  - выходной (отфильтрованный) вектор данных;
///
void CFilter::getConvolution(const std::vector<double>& inData, std::vector<double>& outData)
{
    const int N = inData.size();
    double dx = 1.0/(N-1);
    /// генерируем тестовые данные
//    double xx = 0.0;
//    const double k = 0.68, b = 0.34;
//    for (uint i = 0; i < N; i++)
//    {
//        inData[i] = k*xx + b + 0.1*rnd();
//        xx += dx;
//    }
    /// очищаем очередь
    reset();
    /// центральная часть - фильтруется стандартным (быстро) образом
    /// профильтрованный диапазон данных: outData [shift ... N-1-shift]
    /// вне него данные повторяют исходные
    for (uint i = 0; i < N; i++)
    {
        push( inData[i] );
        outData[i] = inData[i];
        if (i >= shift+shift)
            outData[i-shift+1] = getConvolution();
    }
//    return;
    /// теперь фильтруем левый и правый края полиномом этого же порядка
    std::vector<double> coef(m+1);// вектор коэффициентов полинома
    /// длины нефильтруемых стандартным образом данных с левого и правого
    /// краёв всего массива = shift
    double* x = new double [shift];// псевдоабсциссы
    double* y = new double [shift];
    /// левый край
    for (uint i = 0; i < shift; i++)
    {
        x[i] = (i > 0) ? x[i-1]+dx : -static_cast<int>(shift)*dx;
        y[i] = inData[i];
    }
    coef[0] = outData[shift];
    coef[1] =.5*(-3*outData[shift] + 4*outData[shift+1] - outData[shift+2])/dx;
//    cout << coef[0] << "\t" << coef[1] << endl;
    getCoeffSmoothPoly(coef, x, y, shift);
//    for (uint i = 0; i < coef.size(); i++)
//        cout << coef[i] << "\t";
//    cout << endl;
    /// вычисляем outData [0 ... shift-1]
    for (uint i = 0; i < shift; i++)
        outData[i] = Polynomial( x[i], coef );
    /// правый край
    int iS = N-shift;
    for (uint i = 0; i < shift; i++)
    {
        x[i] = (i > 0) ? x[i-1]+dx : dx;
        y[i] = inData[i+iS];
    }
    coef[0] = outData[iS-1];
    coef[1] =.5*(3*outData[iS-1] - 4*outData[iS-2] + outData[iS-3])/dx;
//    cout << coef[0] << "\t" << coef[1] << endl;
    getCoeffSmoothPoly(coef, x, y, shift);
//    for (uint i = 0; i < coef.size(); i++)
//        cout << coef[i] << "\t";
//    cout << endl;
    /// вычисляем outData [N-shift ... N-1]
    for (uint i = 0; i < shift; i++)
        outData[i+iS] = Polynomial( x[i], coef );
    delete [] x;
    delete [] y;
}
