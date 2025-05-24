#ifndef CCalc_H
#define CCalc_H

#include <QObject>
#include <QVector>
#include <QPointF>
#include <iostream>
//#include "../../Arxiv/cfilter.h"
#include "cfilter.h"

using namespace std;

// глобальные константы:
const double ms = 1e-3;//[s] миллисекунда (ед. времени в данной модели)
///
/// класс: Точка
///
struct CPoint
{
    double x;
    double y;
};
///
/// класс: Граница отображаемой области
///
struct CBounds
{
    double Xmin;
    double Xmax;
    double Ymin;
    double Ymax;
};
///
/// Режимы работы вычислителя:
///     Default - фильтрация зашумленного сигнала заданным фильтром;
///     SyncDetect - синхронное детектирование огибающей высокочастотного АМ сигнала;
///
enum mode {Default, SyncDetect};
///
/// класс моделирующий зашумленный сигнал и его фильтрацию
///
class CCalc : public QObject
{
    Q_OBJECT

public:
    explicit CCalc(QObject *parent = nullptr);
    ~CCalc();
    // 0 - чистый; 1 - зашумленный; 2 - отфильтрованный сигналы
    QVector<QPointF> data[3];// вектора сигналов
    std::vector<CPoint> fixedQueue[3];// очереди сигналов с фиксированной длиной
    std::vector<CPoint> envelope[2];// детектированные сигналы
    unsigned int N;
    unsigned int index;// глобал. индекс отсчётов /INT_MAX = 2147483647/
    unsigned int Nc[2];// верхние значения для массива счётчиков
    unsigned int Nt;// число отсчётов в векторах сигналов
    int idTimer1, idTimer2, iDelay;
    CBounds Bounds;// граница отображаемой области
    // фильтр
    CFilter filter;
    mode SimulationMode = SyncDetect;
    // параметры входного сигнала
    double T, f0, w0;
    std::vector<double> Amp, w, phase;// спектр сигнала
    double SQuadSignal;// сумма квадратов текущей выборки (очереди) сигнала
//    double cDetectedSignal;// детектированная огибающая сигнала (текущее значение)
    unsigned int NT;// число отсчётов в детектируемом сигнале
    double AmpMax;// максим. возможное значение сигнала
    unsigned int Nh;// число гармоник (размерность Amp,...)
    double NoiseAmp;// амплитуда шума
    double SumSTD[3];
    bool isTerminated;
    unsigned int delay;
    unsigned int Nblocks;
    //
    double dt, t, Dt, shift_t;
    // вычисление отсчётов зашумленного и отфильтрованного сигналов
    // и сохранение их в очередь
    void PushSample();
    // приготовить данные для отображения на графике и послать сигнал
    void SendConvertedData();
    // инициализировать процесс моделирования
    void InitializeProcess();

private:
    double VirtualSignalAmp(double cAmp);
    void GenerateSpectr();

signals:
    void UpdatePlotter();
    void ProcessFinished();

public slots:
    void Process();
    void setNoiseAmp(int);
    void setDelay(int);
    void setNblocks(int);
    void setSimulationMode(int);
};

#endif // CCalc_H
