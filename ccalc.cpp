#include "ccalc.h"
#include <cmath>
#include <QDebug>
#include "unistd.h"
#include <QElapsedTimer>

#define RedColor "\033[1;31m"
#define GreenColor "\033[1;32m"
#define YellowColor "\033[1;33m"
#define BlueColor "\033[1;34m"
#define MagentaColor "\033[1;35m"
#define CyanColor "\033[1;36m"
#define ResetColor "\033[0m"

CCalc::CCalc(QObject *parent) :  QObject(parent)
{
    Nblocks = 50;// число обрабатываемых блоков данных между перерисовкой
//    CFilter filt(3, 35, 1, "", false);//13-  25- 27- 29- 31- 33- 35- 37
//    filt.ShowFilterInfo();
    filter.genNewFilter(3, 2*3+1, 3);
    filter.showFilterInfo();
}

//
CCalc::~CCalc()
{
}


//
void CCalc::GenerateSpectr()
{
    // создаём вектора для параметров спектра
    Amp.resize( Nh );
    w.resize( Nh );
    phase.resize( Nh );
    for (unsigned int i = 0; i < Nh; i++)
    {
        Amp[i] = VirtualSignalAmp(1.0);
//        Amp[i] = fabs( 1.0*rand() / RAND_MAX );// 0 <= Amp[i] <= 1.0
        w[i] = 1*(1.0 + 9.0*rand() / RAND_MAX)*w0;// f0 <= f[i] <= 10*f0
        phase[i] = (2.0*rand() / RAND_MAX - 1.0)*M_PI;// -pi <= Amp[i] <= pi
    }
}
///
void CCalc::InitializeProcess()
{
    // задаём параметры спектра сигнала
    f0 = 5.0;//[kHz] основная частота колебаний
    w0 = 2*M_PI*f0;//[rad/s]
    T = 1/f0;//[ms]
    dt = 0.005*T;//[ms] it means freq_discretization = 1/dt = 1/0.005/T = 200*f0 = 1[MHz]; freq_N = 0.5[MHz]
    Nt = static_cast<unsigned int>(3*T/dt);// число будем отображать ~ 3 периода временного отклика
    Dt = Nt*dt;// отображаемый временной интервал
    // генерируем спектр сигнала
    GenerateSpectr();
    if (SimulationMode == SyncDetect)// режим моделирования - синхронное детектирование
    {
        AmpMax = 1.0;// возможное максимальное значение сигнала
        Nh = 1;// гармоника д.б. только одна !
        N = 1000;// кол-во данных в кривых
        for (unsigned int n = 0; n < 2; n++)
            envelope[n].resize( N );
    }
    else// режим моделирования - сглаживание данных (неограниченное во времени)
    {
        Nh = 5;// число гармоник
        // определяем потенциально возможное максимальное значение сигнала
        AmpMax = Amp[0];
        for (unsigned int i = 1; i < Nh; i++)
            AmpMax += Amp[i];
    }
    qDebug()<<"the maximal value of noise ="<<AmpMax*NoiseAmp;
    // устанавливаем верхнее значение для счётчиков counter[0] и counter[1]
    Nc[0] = 2*Nt;// интервал обновления спектра в отсчётах
    Nc[1] = Nblocks*Nt ;// интервал между актами перерисовки кривых на графике в отсчётах
    qDebug()<<"Nc[1] = "<<Nc[1]<<"samples per"<<Nc[1]*dt<<"[virtual ms] are processed between drawing actions\n"
            <<"Nblocks = [Nc[1]/Nt] ="<<Nc[1]/Nt<<"  the number of samples:  Nt ="<<Nt;
    // старт рандомизации
    srand( static_cast<unsigned int>( time(nullptr) ) );
    // создаём фильтр и подготавливаем к работе
    // варианты: {2: 5..11}; {3: 5..11}; {4: 7..11}; {5: 7}
    shift_t = filter.shift*dt;// коррекционный сдвиг фильтрованного сигнала
    filter.reset();// очистка очереди входных данных в фильтре
    qDebug()<<"the time shift between the filtered sample and current ones: shift_t ="<<shift_t;
    // создаём вектора данных и очередей
    for (unsigned int n = 0; n < 3; n++)
    {
        data[n].resize( static_cast<int>(Nt) );
        fixedQueue[n].resize( Nt );
//        qDebug()<<"address data[n] = "<<&data[n];
    }
    // коррекция для правильного отображения начального участка фильтров. сигнала
    for (unsigned int n = 0; n < 3; n++)
        for (unsigned int i = 0; i < Nt; i++)
        {
            if (n == 2)
                fixedQueue[n][i].x =-shift_t;
            else
                fixedQueue[n][i].x = 0;
            fixedQueue[n][i].y = 0;
        }
//    for (index = 0; index < Nt; index++)// тест вычисления лок. индекса - взятия остатка от деления
//        qDebug()<<index<<"  "<<(index % Nt);
//    abort();
}
// вычисляем текущие отсчёты: 0-чистого, 1-зашумленного и 2-отфильтрованного сигналов
// и задвигаем их в очередь /с фиксированной длиной/
void CCalc::PushSample()
{
    // текущие отсчёты: чистый, зашумленный, зашумленный(сдвинутый) и отфильтрованный сигналы
    double oSignal, Signal = 0.0, pSignal, fSignal, noise;
    unsigned int cindex = index % Nt;
    //
    t += dt;
    for (unsigned int i = 0; i < Nh; i++)
        Signal += Amp[i]*sin(w[i]*t + phase[i]);
    Signal += 5.0;
    // чистый сигнал (без шума)
    fixedQueue[0][cindex].x = t;
    if (SimulationMode == SyncDetect)
    {
        fixedQueue[0][cindex].y = Amp[0];// в случае синхронного детектирования
//        AmpMax =  Amp[0];// вариант мультипликативного шума
    }
    else
        fixedQueue[0][cindex].y = Signal;
    noise = NoiseAmp*(2.0*rand() / RAND_MAX - 1.0)*AmpMax;// шум = NoiseAmp[%] -10%..+10%
    Signal += noise;
    // сигнал с шумом
    fixedQueue[1][cindex].x = t;
    fixedQueue[1][cindex].y = Signal;
    // извлекаем последний отсчёт из текущей очереди фильтра зашумленного сигнала
    pSignal = filter.GetLastSampleFromQueue();
    // фильтрация
    filter.push( Signal );
    fSignal = filter.getConvolution();
    // реализован оптимизированный вариант расчёта квадрата локальной выборки:
    // удаляем квадрат уходящего отсчёта из суммы квадратов текущей выборки сигнала
    if (SimulationMode == SyncDetect)
        SQuadSignal -= 2*fixedQueue[2][cindex].y*fixedQueue[2][cindex].y;
    // отфильтрованный сигнал
    fixedQueue[2][cindex].x = t - shift_t;
    fixedQueue[2][cindex].y = fSignal;
    // добавляем квадрат пришедшего отсчёта из суммы квадратов текущей выборки сигнала
    if (SimulationMode == SyncDetect)
        SQuadSignal += 2*fixedQueue[2][cindex].y*fixedQueue[2][cindex].y;
    // стандартное отклонение зашумлённого от исходного сигнала
    SumSTD[0] += noise*noise;
    // стандартное отклонение отфильтрованного от зашумлённого сигнала
    SumSTD[1] += (fSignal - pSignal)*(fSignal - pSignal);
    // стандартное отклонение отфильтрованного от исходного сигнала
    // характеризует точность восстановления исходного сигнала
    if (cindex >= filter.shift)
        oSignal = fixedQueue[0][cindex-filter.shift].y;
    else
        oSignal = fixedQueue[0][Nt-(filter.shift-cindex)].y;
    SumSTD[2] += (fSignal - oSignal)*(fSignal - oSignal);
//    if (index > 1000 && cindex > filter.shift+5)
//    {
//        qDebug()<<oSignal<<"  "<<Signal<<"  "<<fSignal;
//        for (uint i = cindex; i > cindex-filter.shift; i--)
//            qDebug()<<i<<":"<<fixedQueue[0][i].y<<"  "<<fixedQueue[1][i].y<<"  "<<fixedQueue[2][i].y;
//        qDebug()<<"++++++++++++++++++++++++++++++++++++++++";
//    }
    // определяем текущие границы графиков
    Bounds.Xmin = t - Dt;
    Bounds.Xmax = t;
    if (Bounds.Ymin > Signal)
        Bounds.Ymin = Signal;
    if (Bounds.Ymax < Signal)
        Bounds.Ymax = Signal;
}
//
void CCalc::Process()
{
    QElapsedTimer timer;
    long long StartTime, WorkTime;
    unsigned int counter[2] = {0, 0};// массив счётчиков
    unsigned int j = 0;// индекс массива огибающих

    if (SimulationMode == SyncDetect) // режим синхронного детектирования
    {
        NT = static_cast<unsigned int>( round( 2.0*2*M_PI/w[0]/dt ) );//д.б. после ResetSimulation()
    //    qDebug()<<"T="<<T<<" 2*M_PI/w[0]="<<2*M_PI/w[0]<<"  dt="<<dt<<"  NT="<<NT;
    }
    else // режим сглаживания
    {
        N = 1;// для выполнения не актуального в этом режиме условия j < N
    }
    isTerminated = false;
    // задаем начальные условия
    SumSTD[0] = SumSTD[1] = SumSTD[2] = 0.0;
    index = 0;
    t =-dt;
    Bounds.Ymin = Bounds.Ymax = 0;//мин. и макс. значения динамич. переменной на графике
    qDebug()<<RedColor<<"\n  The process is starting..."<<ResetColor;
    timer.start();
    StartTime = timer.elapsed();
    while (j < N && !isTerminated)
    {
        // загрузка отсчётов сигнала (исходного/зашумленного/отфильтрованного) в очереди
        PushSample();
        // синхронное детектирование огибающей сигнала & изменение спектра сигнала "на лету"
        if (counter[0] < Nc[0])
            counter[0]++;
        else
        {
            counter[0] = 0;
            // синхронное детектирование огибающей сигнала
            if (SimulationMode == SyncDetect)
            {
                // сохранение исходной и детектированной огибающих
                envelope[0][j].x = t;
                envelope[0][j].y = Amp[0];
                envelope[1][j].x = t - shift_t;
                envelope[1][j].y = sqrt( SQuadSignal / NT );//синхронно детектированный сигнал (огибающая)
                j++;
            }
            // трансформация спектра сигнала
            for (unsigned int i = 0; i < Nh; i++)
            {
                Amp[i] = VirtualSignalAmp(Amp[i]);
                w[i] += 1e-9*rand()/RAND_MAX*w0;
                phase[i] = 1e-2*(2.0*rand() / RAND_MAX - 1.0)*M_PI;
            }
        }
        // динамическое отображение данных в плоттере
        if (counter[1] < Nc[1]-1)
            counter[1]++;
        else // конвертация очереди в вектор и отображение на графике
        {
            counter[1] = 0;
            SendConvertedData();
            // отображение STD сигнала - "оценки" шума
            qDebug()<<index<<"STD : {Sig-oSig} ="<<sqrt( SumSTD[0] / index )
                    <<" {fSig-Sig} ="<<sqrt( SumSTD[1] / index )
                    <<" {fSig-oSig} /recovery accuracy/ ="<<sqrt( SumSTD[2] / index );
            usleep( delay );//delay[msec]
        }
        index++;
    }
    WorkTime = timer.elapsed() - StartTime;//[ms]
    double Nsamples = index/1e6;//[Msamples]
    unsigned int cNblocks = static_cast<unsigned int>(index)/Nc[1];
    double SleepTime = cNblocks*delay*1e-3;//[ms]
    qDebug()<<CyanColor<<"\n  The process is stopped !\n"<<ResetColor;
    qDebug()<<"the "<<index<<"samples were received for"<<WorkTime*ms<<"[s] ("<<WorkTime<<"[ms] )";
    qDebug()<<cNblocks<<"blocks were processed";
    qDebug()<<"the effective rate is a"<<Nsamples/(WorkTime*ms)<<"[Msamples/s]";
    qDebug()<<"the time of calcul-s & drawing plots is ~"<<(WorkTime-SleepTime)<<"[ms]";
    qDebug()<<"the sleeping time is ~"<<SleepTime<<"[ms]";
    qDebug()<<"the rate (without sleep time) is a"<<Nsamples/((WorkTime-SleepTime)*ms)<<"[Msamples/s]";
    if (SimulationMode == SyncDetect) // синхронное детектирование
    {
        if (j == N)
            emit ProcessFinished();
    }
}

// вычисление мгновенной амплитуды данной гармоники виртуального сигнала
// cAmp - текущее значение амплитуды этой гармоники
double CCalc::VirtualSignalAmp(double cAmp)
{
    double Amp;
    if (SimulationMode == SyncDetect)
    {
        double t = index*dt, tend = Nc[0]*N*dt;
        double ct = t-tend/2;// текущий момент времени отн-но середины интервала
        double cnt = ct/(0.15*tend);// текущий нормированный момент времени
        Amp = exp(-cnt*cnt) + 0.05*(cos(2*M_PI/(tend/20)*ct) + 1.0);
    }
    else // 0 <= Amp[i] <= 1.0
    {
        double dAmp = 1e-3*(2.0*rand() / RAND_MAX - 1.0 );// малое случайное приращение амплитуды
        if (cAmp + dAmp < 1.0)  Amp = cAmp + dAmp;  else Amp = 1.0 - dAmp;
        Amp = fabs( Amp );
    }
    return Amp;
}

// отображаем сигналы на графике
void CCalc::SendConvertedData()
{
    // конвертация данных из очереди в нормальный вектор
    if (Nc[1] % Nt == 0) // I вариант - для Nc[1] кратного Nt
    {
        for (unsigned int n = 0; n < 3; n++)
        {
            for (unsigned int i = 0; i < Nt; i++)
                data[n].replace( static_cast<int>(i), QPointF(fixedQueue[n][i].x, fixedQueue[n][i].y) );
        }
    }
    else // II вариант - для произвольного Nc[1] не кратного Nt ???
    {
        unsigned int ii = index % Nt;
        for (unsigned int n = 0; n < 3; n++)
        {
            for (unsigned int i = ii; i < Nt; i++)//ii..Nt-1 = 0..Nt-ii-1
                data[n].replace( static_cast<int>(i-ii), QPointF(fixedQueue[n][i].x, fixedQueue[n][i].y) );
            for (unsigned int i = 0; i < ii; i++)//0..ii-1 = Nt-ii..Nt-1
                data[n].replace( static_cast<int>(i-ii+Nt), QPointF(fixedQueue[n][i].x, fixedQueue[n][i].y) );
        }
    }
//    double STD[3] = {0, 0, 0};
//    for (unsigned int i = 0; i < Nt; i++)
//    {
//        STD[0] += pow(fixedQueue[1][i].y - fixedQueue[0][i].y, 2);//зашум. - исход.
//        if (i < Nt-filter.shift)
//        {
//            unsigned int ii = i+filter.shift;
//            STD[1] += pow(fixedQueue[2][i].y - fixedQueue[0][ii].y, 2);//фильтров. - исход.
//            STD[2] += pow(fixedQueue[2][i].y - fixedQueue[1][ii].y, 2);//фильтров. - зашум.
////            qDebug()<<fixedQueue[0][i].y<<"  "<<fixedQueue[1][i].y<<"  "<<fixedQueue[2][ii].y;
//        }
//    }
//    qDebug()<<STD[2]/STD[0];//<<" "<<STD[1]<<" "<<STD[2];
//    qDebug()<<Bounds.Xmin<<"  "<<Bounds.Xmax<<"  "<<Bounds.Ymin<<"  "<<Bounds.Ymax;

    emit UpdatePlotter();
}

// установка амплитуды шума
void CCalc::setNoiseAmp(int _NoiseAmp)
{
    NoiseAmp = 0.001*_NoiseAmp;
//    qDebug()<<"NoiseAmp = "<<NoiseAmp;
}

//
void CCalc::setDelay(int Delay/*[ms]*/)
{
    delay = static_cast<unsigned int>(Delay)*1000;//[us]
//    qDebug()<<"delay = "<<delay/1000;
}

//
void CCalc::setNblocks(int _Nblocks)
{
    // тонкая подстройка: для корректного рисования
    Nblocks = static_cast<unsigned int>(_Nblocks);
}

// установка режима работы
void CCalc::setSimulationMode(int _SimulationMode)
{
    SimulationMode = static_cast<mode>(_SimulationMode);
    switch (_SimulationMode)
    {
    case SyncDetect:
        Nh = 1;
        qDebug()<<"the simulation mode: sync detect";
        break;
    default:// Default
        Nh = 5;
        qDebug()<<"the simulation mode: default - smoothing";
        break;
    }
}
