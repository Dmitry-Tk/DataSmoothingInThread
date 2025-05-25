/*! \file projects/Arxiv/cfilter.h */
#ifndef CFILTER_H
#define CFILTER_H
#include <vector>
#include <QString>

/**
 * !\class CFilter
 *
 * Класс - многокаскадный фильтр позволяет сглаживать (фильтровать) поток данных
 * как "на лету", так и пакетами. Заполнение фильтра данными происходит конвейерным образом,
 * т.е. без перезаписи содержимого.\n
 * Однокаскадный фильтр:\n
 *     newSample -> in[0]->in[1]...->in[length-2]->in[length-1]-x\n
 * Многокаскадный фильтр:\n
 *  0 уровень:*           newSample*-> In[0][0]->In[0][1] *...* ->In[0][length-2]->In[0][length-1]-x\n
 *  1 уровень:*  _getConvolution(0)*-> In[1][0]->In[1][1] *...* ->In[1][length-2]->In[1][length-1]-x\n
 *      ...   *            ...     *                      *...* \n
 *  N уровень:*_getConvolution(N-1)*-> In[N][0]->In[N][1] *...* ->In[N][length-2]->In[N][length-1]-x\n
 *          N = level-1 - индекс последнего каскада.\n
 *
 * результат фильтрации: getConvolution()
**/

class CFilter
{
public:
    /// \brief конструктор
    CFilter(unsigned int m = 3, unsigned int length = 7, unsigned int N = 1,
            QString FilterName = "");
    /// \brief деструктор
    ~CFilter();
    /// \brief генератор нового фильтра
    void genNewFilter(unsigned int m = 3 , unsigned int length = 7, unsigned int N = 1, QString FilterName = "");
    /// \brief ввод нового отсчёта в очередь фильтра
    void push(double newSample);
    /// \brief очистка очереди (входного потока данных) фильтра
    void reset();
    /// \brief вычисление текущего фильтруемого отсчёта, т.е. вычисление свёртки фильтра и текущих входных данных
    double getConvolution();
    /// \brief выполнение фильтрации входного вектора данных, т.е. вычисление его свёртки с фильтром
//    void getConvolution(const std::vector<double>& inData, std::vector<double>& outData);
    void getConvolution(const std::vector<double>& inData, std::vector<double>& outData);
    /// \brief информация о фильтре
    void showFilterInfo(bool in_detail = true);
    /// \brief получить параметры фильтра
    void getFilterParams(unsigned int& m, unsigned int& length, unsigned int& N);
    /// \brief возвращение текущего последнего отсчёта из очереди фильтра
    double GetLastSampleFromQueue();

    unsigned int m;///< кол-во коэффициентов в аппроксимационном полиноме
    unsigned int halflength;///< половина длины фильтра
    unsigned int shift;///< сдвиг между входным и фильтруемым отсчётом
    unsigned int N;///< число каскадов фильтра
    static unsigned int Nfilters;///< число фильтров

private:
    bool ShowObjectInfo;///< ключ: определяет отображать или нет состояние объекта в процессе работы (для отладки)
    int length;///< длина фильтра - число точек, по которым вычисляется отфильтрованное значение
    /// \brief глобальный индекс нового добавленного в фильтр отсчёта
    /// \details указывает сколько отсчётов уже прошло через фильтр
    long long index;
    /// \brief локальный индекс нового добавленного в фильтр отсчёта
    /// \details указывает индекс последнего добавленного отсчёта в очереди фильтра
    int cindex;
    /// \brief используется и в одно-, и многокаскадном фильтре
    /// (в последнем случае для сопоставления текущему отфильтрованному отсчёту его нефильтрованного значения)
    /// \details однокаскадная очередь:
    ///     \f$in_{0}\ -\ last\ sample\quad\cdots\quad in_{length-1}\ -\ first\ sample\f$
    double* in = nullptr;
    /// \brief используется только в многокаскадном фильтре
    /// \details многокаскадная очередь:
    ///     - 1 каскад: \f$\{In_{0,0}\ -\ last\ sample\quad\cdots\quad In_{0,length-1}\ -\ first\ sample\}\f$
    ///     - 2 каскад: \f$\{In_{1,0}\ -\ last\ sample\quad\cdots\quad In_{1,length-1}\ -\ first\ sample\}\f$
    ///     -               ...
    ///     - N каскад: \f$\{In_{N-1,0}\ -\ last\ sample\quad\cdots\quad In_{N-1,length-1}\ -\ first\ sample\}\f$
    double** In = nullptr;
    double* coeff = nullptr;///< коэффициенты фильтра
    QString FilterName;///< название фильтра
    double thSTD;///< теорет. оценка STD при фильтрации
    /// \brief вычисление коэффициентов фильтра
    void getFLoatCoeff();
    /// \brief вычисление свёртки n-го каскада фильтра
    double _getConvolution(unsigned int level);
};

#endif // CFILTER_H
