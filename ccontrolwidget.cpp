#include "ccontrolwidget.h"
#include <QString>
#include <QtCharts/QLogValueAxis>
#include <QtCharts/QValueAxis>
#include <QFileDialog>
#include <QVector>
//#include <../../Arxiv/mathlib.h>
#include "mathlib.h"


using namespace QtCharts;

static const QColor colorID[12] = { Qt::red, Qt::blue, Qt::green, Qt::magenta, Qt::cyan,
                                    Qt::yellow, Qt::darkGreen, Qt::darkRed, Qt::darkBlue,
                                    Qt::darkMagenta, Qt::darkCyan, Qt::darkYellow };

CControlWidget::CControlWidget(QWidget *parent) :  QWidget(parent)
{
    QGroupBox *groupBox1 = new QGroupBox("mode");
    QVBoxLayout *vbox1 = new QVBoxLayout;
    radioBtns.resize(3);
    for (int i=0; i<radioBtns.size(); i++)
    {
        radioBtns[i] = new QRadioButton;
        vbox1->addWidget(radioBtns[i]);
    }
    radioBtns[0]->setText("default");
    radioBtns[1]->setText("sync detect");
    radioBtns[2]->setText("prototype");// для будущего использования
    radioBtns[2]->setEnabled( false );// блокируем
    groupBox1->setLayout(vbox1);
    ///
    QGroupBox *groupBox2 = new QGroupBox("model parameters");
    QGridLayout *GridLayout = new QGridLayout;
    spinboxNblocks = new QSpinBox;
    QLabel* lab1 = new QLabel("number of blocks");
    spinboxDelay = new QSpinBox;
    QLabel* lab2 = new QLabel("delay after drawing");
    checkboxYScale = new QCheckBox("Y-axis: log scale/normal");
    GridLayout->addWidget(lab1,0,0); GridLayout->addWidget(spinboxNblocks,0,1);
    GridLayout->addWidget(lab2,1,0); GridLayout->addWidget(spinboxDelay,1,1);
    GridLayout->addWidget(checkboxYScale,2,0);
    groupBox2->setLayout(GridLayout);
    spinboxDelay->setRange(1, 1000);
    spinboxDelay->setValue(100);
    spinboxDelay->setSuffix( tr(" [ms]") );
    spinboxDelay->setAlignment(Qt::AlignCenter);
    spinboxNblocks->setRange(50, 1000);
    spinboxNblocks->setValue(40);
    spinboxNblocks->setAlignment(Qt::AlignCenter);
    checkboxYScale->setCheckable( true );
    checkboxYScale->setChecked( false );
    ///
    //    QGroupBox *groupBox3 = new QGroupBox("controls");
    QVBoxLayout *vbox3 = new QVBoxLayout;
    btnStartStop = new QPushButton("Start&&Pause Calc-s");
    btnSmoothData = new QPushButton("Smooth Data");
    btnResetChart = new QPushButton("Reset Chart");
    btnSetOriginalZoom = new QPushButton("Original Zoom");
    btnSetFilterParams = new QPushButton("Filter param-s");
    vbox3->addWidget(btnStartStop);
    vbox3->addWidget(btnSetFilterParams);
    vbox3->addWidget(btnSmoothData);
    vbox3->addWidget(btnSetOriginalZoom);
    vbox3->addWidget(btnResetChart);    
    //    groupBox3->setLayout(vbox3);
    btnStartStop->setFlat(true);
    btnStartStop->setCheckable(true);
    btnSmoothData->setFlat(true);
    btnSetOriginalZoom->setFlat(true);
    btnResetChart->setFlat(true);
    btnSetFilterParams->setFlat(true);
//    btnSetFilterParams->setFlat(true);

    sliderNoise = new QSlider;
    sliderNoise->setOrientation(Qt::Vertical);
    sliderNoise->setMinimum(0);
    sliderNoise->setMaximum(1000);
    sliderNoise->setValue(0);

    QHBoxLayout *hbox = new QHBoxLayout;
    hbox->addWidget(sliderNoise);
    hbox->addWidget(groupBox1);
    hbox->addStretch(1);
    hbox->addWidget(groupBox2);
    hbox->addStretch(1);
    hbox->addLayout(vbox3);
    //    hbox->addWidget(groupBox3);
    ///
    setLayout(hbox);
    ///
    InitPlotter();
    SmoothingObject.moveToThread( &CalcThread );
    ///
    /// соединение сигналов и слотов
    /// старт - стоп вычислений
    connect( btnStartStop, SIGNAL( clicked() ), this, SLOT( StartStopCalc() ) );
    connect( &CalcThread, SIGNAL( started() ), &SmoothingObject, SLOT( Process() ) );
    connect( &SmoothingObject, SIGNAL( UpdatePlotter() ), this, SLOT( onUpdatePlotter() ) );
    connect( &SmoothingObject, SIGNAL( ProcessFinished() ), this, SLOT( StartStopCalc() ) );
    connect( &SmoothingObject, SIGNAL( ProcessFinished() ), this, SLOT( onProcessFinished() ) );
    /// загружаем и сглаживаем загружаемые данные
    connect( btnSmoothData, SIGNAL( clicked() ), this, SLOT( SmoothData() ) );
    /// Y-axis: log/normal
    connect( checkboxYScale, SIGNAL(clicked(bool)), this, SLOT(setYScale(bool)));
    /// очистка плоттера
    connect( btnResetChart, SIGNAL( clicked() ), this, SLOT( ClearChart() ) );
    /// сброс плоттера в исходный масштаб
    connect( btnSetOriginalZoom, SIGNAL( clicked() ), this, SLOT( SetOriginalZoom() ) );
    /// установка NoiseAmp - амплитуды шума
    connect( sliderNoise, SIGNAL( valueChanged(int) ), &SmoothingObject, SLOT( setNoiseAmp(int) ), Qt::DirectConnection );
    connect( sliderNoise, SIGNAL( valueChanged(int) ), this, SLOT( setNoiseAmp(int) ) );
    /// установка Delay - задержки после прорисовки графика
    //    connect( spinboxDelay, &QSpinBox::valueChanged, &SmoothingObject, &CCalc::setDelay);
    connect( spinboxDelay, SIGNAL( valueChanged(int) ), &SmoothingObject, SLOT( setDelay(int) ), Qt::DirectConnection );
    /// установка Nblocks - обрабатываемое число блоков отсчётов перед прорисовкой графика
    connect( spinboxNblocks, SIGNAL( valueChanged(int) ), &SmoothingObject, SLOT( setNblocks(int) ), Qt::DirectConnection );
    /// установка Nblocks - обрабатываемое число блоков отсчётов перед прорисовкой графика
    connect( btnSetFilterParams, SIGNAL( clicked() ), this, SLOT( setFilterParams() ) );

    /// установка режима вычислений: default - smoothing; sync detect - синхронное детектирование
    for (int i=0; i<radioBtns.size(); i++)
        connect( radioBtns[i], SIGNAL(clicked()), this, SLOT(setSimulationMode()));
    connect( this, SIGNAL(SimulationModeChanged(int)), &SmoothingObject, SLOT(setSimulationMode(int)), Qt::DirectConnection );
    ///
    /// начальная инициализация и предустановка параметров
    radioBtns[0]->setChecked(true);    setSimulationMode();
    SmoothingObject.setNoiseAmp( sliderNoise->value() );
    SmoothingObject.setDelay( spinboxDelay->value() );
    SmoothingObject.setNblocks( spinboxNblocks->value() );
    sliderNoise->setValue( 100 );
}
///
CControlWidget::~CControlWidget()
{
    SmoothingObject.isTerminated = true;
    CalcThread.quit();
    CalcThread.wait();
}
///
void CControlWidget::InitPlotter()
{
    ChartView = new QChartView;
    ChartView->resize(800, 600);
    ChartView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    ChartView->setRubberBand(QChartView::RectangleRubberBand);
    chart = new QChart;
    chart->setTheme(QChart::ChartThemeDark);
    ChartView->setChart(chart);
    for (unsigned int n = 0; n < 3; n++)
        series[n] = nullptr;
}
//
void CControlWidget::StartStopCalc()
{
    if (!CalcThread.isRunning())
    {
        btnStartStop->setChecked( true );
        // перед соединением каналов данных вызов
        // SmoothingObject.InitializeProcess() обязателен !!!
        SmoothingObject.InitializeProcess();
        ConnectDataChannel();
        CreateSeries();
        CalcThread.start(QThread::LowPriority);
    }
    else
    {
        btnStartStop->setChecked( false );
        SmoothingObject.isTerminated = true;
        CalcThread.quit();
        CalcThread.wait();
    }
}
//
void CControlWidget::ConnectDataChannel()
{
    Bounds = &(SmoothingObject.Bounds);
    for (int n = 0; n < 3; n++)
    {
        data[n] = &(SmoothingObject.data[n]);
//        qDebug()<<"ConnectDataChannel(): address data[n] = "<<data[n];
    }
}
//
void CControlWidget::CreateSeries()
{
    if (series[0] != nullptr)
        return;
    for (unsigned int n = 0; n < 3; n++)
    {
        series[n] = new QLineSeries();
        series[n]->setColor( colorID[n+3] );
        series[n] = dynamic_cast<QXYSeries*>(series[n]);
        chart->addSeries( series[n] );
    }
    series[0]->setName( tr("original signal") );
    series[1]->setName( tr("noisу signal") );
    series[2]->setName( tr("filtered signal") );
    chart->createDefaultAxes();
    chart->axisX()->setTitleText("t, [ms]");
    chart->axisY()->setTitleText("Signals");
    chart->legend()->setAlignment( Qt::AlignRight );
    for (int i = 0; i < 3; i++)
    {
        series[i]->setVisible( true );
        series[i]->setPen( QPen(colorID[i], 3.0) );
    }
}
// установка режима вычислений
void CControlWidget::setSimulationMode()
{
    for (int n = 0; n < radioBtns.size(); n++)
        if (radioBtns[n]->isChecked() && radioBtns[n]->isEnabled())
            emit SimulationModeChanged( n );
}

// установка уровня шума
void CControlWidget::setNoiseAmp(int NoiseAmp)
{
    sliderNoise->setToolTip( QString::number(0.1*NoiseAmp) + "%" );
}

///
void CControlWidget::setFilterParams()
{
    DefineFilterParamDialog dialog(this);
    if (dialog.exec() == QDialog::Rejected)
        return;
    SFilterParams* FilterParams = dialog.getFilterParams();
    // создаём фильтр и подготавливаем к работе
    // варианты: {2: 5..11}; {3: 5..11}; {4: 7..11}; {5: 7}
    SmoothingObject.filter.genNewFilter(FilterParams->n_poly, 2*FilterParams->halfwidth+1, FilterParams->cascades);
    SmoothingObject.filter.showFilterInfo();
    SmoothingObject.shift_t = SmoothingObject.filter.shift*SmoothingObject.dt;// коррекционный сдвиг фильтрованного сигнала
    SmoothingObject.filter.reset();// очистка очереди входных данных в фильтре
}
//
void CControlWidget::ClearChart()
{
    chart->removeAllSeries();
    for (unsigned int n = 0; n < 3; n++)
        series[n] = nullptr;
    cout<<"All series have been removed !"<<endl;
}

// сброс масштаба плоттера в исходный
void CControlWidget::SetOriginalZoom()
{
    chart->zoomReset();
}

// отображаем сигналы на графике
void CControlWidget::onUpdatePlotter()
{
    for (unsigned int n = 0; n < 3; n++)
        series[n]->replace( *data[n] );

    chart->axisX()->setMin( Bounds->Xmin );
    chart->axisX()->setMax( Bounds->Xmax );
    chart->axisY()->setMin( Bounds->Ymin );
    chart->axisY()->setMax( Bounds->Ymax );
}

// после завершения процесса загружаем исходные и измеренные данные
void CControlWidget::onProcessFinished()
{
    int N = static_cast<int>(SmoothingObject.envelope[0].size());
    QVector<QPointF> EnvData( static_cast<int>(N) );
    Bounds->Xmin = SmoothingObject.envelope[0][0].x;
    Bounds->Xmax = SmoothingObject.envelope[0][static_cast<uint>(N-1)].x;
    Bounds->Ymin = Bounds->Ymax = SmoothingObject.envelope[0][0].y;
    double x, y;
    for (unsigned int n = 0; n < 2; n++)
    {
        for (int i = 0; i < N; i++)
        {
            uint ii = static_cast<uint>(i);
            if (LogScale)
            {
                x = SmoothingObject.envelope[n][ii].x;
                y = log10(SmoothingObject.envelope[n][ii].y);
            }
            else
            {
                x = SmoothingObject.envelope[n][ii].x;
                y = SmoothingObject.envelope[n][ii].y;

            }
            EnvData.replace( i, QPointF(x, y) );
            if (Bounds->Ymin > SmoothingObject.envelope[n][ii].y) Bounds->Ymin = SmoothingObject.envelope[n][ii].y;
            if (Bounds->Ymax < SmoothingObject.envelope[n][ii].y) Bounds->Ymax = SmoothingObject.envelope[n][ii].y;
        }
        series[n]->replace( EnvData );
    }
    series[2]->setVisible( false );
    //    for (unsigned int i = 0; i < static_cast<unsigned int>(N); i++)
    //        qDebug()<<fabs(100*(envelope[1][i].y - envelope[0][i].y)/envelope[0][i].y)<<"%";
    if (LogScale)
    {
        Bounds->Ymin = -4; Bounds->Ymax = 1;
    }
    chart->axisX()->setRange(Bounds->Xmin, Bounds->Xmax);
    chart->axisY()->setRange(Bounds->Ymin, Bounds->Ymax);
    series[0]->setName( tr("the origin envelope") );
    series[1]->setName( tr("the detected envelope") );
    cout<<"SeriesCount = "<<chart->series().size()<<endl;
}

// масштаб Y оси: LogScale = true/false = логарифмический/нормальный
void CControlWidget::setYScale(bool NewLogScale)
{
    LogScale = NewLogScale;
    if (SmoothingObject.SimulationMode == SyncDetect)
    {
        int N = static_cast<int>(SmoothingObject.envelope[0].size());
        QVector<QPointF> data(static_cast<int>(N));
        double x, y, minY, maxY;
        minY = maxY = 0.5;
        for (int n = 0; n < 2; n++)
        {
            for (int k = 0; k < N; k++)
            {
                x = series[n]->points()[k].rx();
                if (LogScale) y = log10( series[n]->points()[k].ry() );
                else y = exp( log(10.0)*series[n]->points()[k].ry() );
                data[k] = QPointF(x, y);
                if (minY > y) minY = y;
                if (maxY < y) maxY = y;
            }
            series[n]->replace(data);
        }
        QValueAxis* YAxis = dynamic_cast<QValueAxis*>(chart->axisY());
        if (LogScale)
        {
            minY = -4; maxY = 1;
            YAxis->setTickCount(6);
        }
        else
            YAxis->setTickCount(5);
        chart->axisY()->setRange(minY, maxY);
    }
    /*
    if (LogScale && SimulationMode == SyncDetect)
    {
        QLogValueAxis* YLogAxis = new QLogValueAxis;
        YLogAxis->setRange(-4.0, 1.0);
        YLogAxis->setMinorTickCount(1.0);
        chart->removeAxis(chart->axisX());
        for (int n = 0; n < 2; n++)
            chart->setAxisY(YLogAxis, series[0]);
    }
    */
}

/// сглаживаем данные фильтром
void CControlWidget::SmoothData()
{
    QString FullFileNameString = QFileDialog::getOpenFileName(this,
            tr("Open Data File"), "/home/Work/Scilab_work/Optical_activity", tr("Data Files (*.dat *.txt)"));
    if (FullFileNameString == "")
        return;
    cout<<"the "<<FullFileNameString.toLatin1().data()<<" file\n\thas been loaded !"<<endl;
    QFileInfo FileInfo = QString( FullFileNameString );
    QString FileName = QString( FileInfo.fileName() );
    std::vector<QPointF> Data;
    LoadData( FullFileNameString, Data );
    if (Data.size() == 0)
    {
        cout<<"the data file is empty !"<<endl;
        return;
    }
    for (uint i = 0; i < Data.size(); i++)
        Data[i].setX( static_cast<qreal>(i) );
    /// Сглаживание данных
    std::vector<QPointF> SmoothedData;
    /// I.
    SmoothSeriesI(Data, SmoothedData);
    /// II.
//    SmoothSeriesII(Data, SmoothedData);
    /// III.
//    SmoothSeriesIII(Data, SmoothedData);
    /// загружаем исходные данные в график
//    QLineSeries* origSeries = new QLineSeries;
    QScatterSeries* origSeries = new QScatterSeries;
    origSeries->replace( QVector<QPointF>::fromStdVector( Data ) );
    origSeries->setName( FileName );
    origSeries->setMarkerSize( 1.0 );
    chart->addSeries( origSeries );
    /// загружаем сглаженные данные в график
    if (SmoothedData.size() == 0)
    {
        cout<<"the data smoothing has been unsuccessfull !"<<endl;
        return;
    }
//    for (uint i = 0; i < SmoothedData.size(); i++)
//        cout<<i<<" : x = "<<SmoothedData[i].x()<<"; y = "<<SmoothedData[i].y()<<endl;
    /// загружаем отфильтрованные данные в график
    QLineSeries* smoothedSeries = new QLineSeries;
    smoothedSeries->replace( QVector<QPointF>::fromStdVector( SmoothedData ) );
    smoothedSeries->setName(  tr("Smoothed ") + FileName );
    chart->addSeries( smoothedSeries );
    chart->createDefaultAxes();

    // emit DataLoaded( static_cast<QXYSeries*>(series) );
}

/// I. вариант сглаживания данных:
/// прогонка данных через полиномиальный многокаскадный фильтр
void CControlWidget::SmoothSeriesI(std::vector<QPointF>& Data, std::vector<QPointF>& SmoothedData)
{
    /// создаём фильтр и подготавливаем к работе
    /// варианты: {2: 5..11}; {3: 5..11}; {4: 7..11}; {5: 7}
    CFilter filter;
    filter.genNewFilter(3, 2*80+1, 1);
    std::vector<double> Ydata(Data.size()), filtYdata(Data.size());
    for (uint i = 0; i < Data.size(); i++)
        Ydata[i] = Data[i].y();
    /// фильтруем данные
    filter.getConvolution(Ydata, filtYdata);
    double shift_lambda = filter.shift*(Data[1].x() - Data[0].x());
    SmoothedData.resize( Data.size() );
    for (uint i = 0; i < Data.size(); i++)
    {
        Data[i] = QPointF( Data[i].x() - 0*shift_lambda, Ydata[i] );
        SmoothedData[i] = QPointF( Data[i].x() - 0*shift_lambda, filtYdata[i] );
    }
}

/// II. вариант сглаживания данных:
/// простое усреднение по ближайшим значениям
void CControlWidget::SmoothSeriesII(const std::vector<QPointF>& Data, std::vector<QPointF>& SmoothedData)
{
    unsigned long N = Data.size();
    uint startIndex, stopIndex, Np = 60, Nb = Np/3;//Nb <= Np/2
    for (uint i = 0; i < N; i++)
        if (i % Np == 0)
        {
            if (i > Nb-1)  startIndex = i - Nb; else startIndex = 0;
            if (i < N-Nb)  stopIndex = i + Nb;  else stopIndex = static_cast<uint>(N)-1;
            double meanY = 0;
            for (uint j = startIndex; j < stopIndex+1; j++)
                meanY += Data[j].y();
            meanY /= stopIndex - startIndex + 1;
            SmoothedData.push_back( QPointF(Data[i].x(), meanY) );
        }
}

/// III. вариант сглаживания данных:
/// полиномиальное сглаживание старая версия, нужна ли ?! не нужна
void CControlWidget::SmoothSeriesIII(const std::vector<QPointF>& Data, std::vector<QPointF>& SmoothedData)
{
    uint N = static_cast<uint>(Data.size()), Ns = N;
    double *xs, *ys, *x, *y;
    double Xmin = Data[0].x(), Xmax = Data[N-1].x();
    bool XDirection;

    xs = new double[Ns];
    ys = new double[Ns];
    x = new double[N];
    y = new double[N];
    if (Xmin < Xmax)
        XDirection = true;//straight direction X
    else
        XDirection = false;//reverse direction X
    if (!XDirection)
    {
        double p = Xmax;
        Xmax = Xmin;
        Xmin = p;
    }
    for (unsigned int i = 0; i < N; i++)
    {
        if (XDirection)
        {
            x[i] = (2*Data[i].x() - (Xmax - Xmin))/(Xmax - Xmin);
            y[i] = Data[i].y();
        }
        else
        {
            x[N-1-i] = (2*Data[i].x() - (Xmax - Xmin))/(Xmax - Xmin);
            y[N-1-i] = Data[i].y();
        }
    }
//    PolynomSmoothing(xs, ys, Ns, x, y, N, 10);
    InterpCompression(xs, ys, Ns, 15, x, y, N);
    SmoothedData.resize( N );
    for (unsigned int i=0; i < N; i++)
        SmoothedData[i] = QPointF((Xmax - Xmin)/2 + xs[i]*(Xmax - Xmin)/2, ys[i]);

    delete[] xs;
    delete[] ys;
    delete[] x;
    delete[] y;
}

/// загружаем данные из файла
void CControlWidget::LoadData(const QString& FullFileNameString, std::vector<QPointF>& data)
{
    QFile file( FullFileNameString );
    if (file.open(QIODevice::ReadOnly))
    {
        QTextStream in( &file );
        qreal x, y;
        in>>x>>y;
        do
        {
            data.push_back( QPointF(x, y) );
            in>>x>>y;
        }
        while (!in.atEnd());
        file.close();
    }
}

/// сохраняем данные в файл
void CControlWidget::SaveData(const QString& FullFileNameString, const std::vector<QPointF>& Data)
{
    QFile file( FullFileNameString );
    std::vector<QPointF> data;
    if (file.open(QIODevice::WriteOnly))
    {
        QTextStream out( &file );
        for (uint i = 0; i < Data.size(); i++)
            out<<Data[i].x()<<"\t"<<Data[i].y();
        file.close();
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \brief The DefineAbsorpParamDialog class
///
/// рабочие варианты: {2: 5..11}; {3: 5..11}; {4: 7..11}; {5: 7}
int DefineFilterParamDialog::n_poly = 3;
int DefineFilterParamDialog::halfwidth = 3;
int DefineFilterParamDialog::cascades = 3;
DefineFilterParamDialog::DefineFilterParamDialog(QWidget *parent)
    : QDialog(parent),
    m_n_poly( new QSpinBox ),
    m_halfwidth( new QSpinBox ),
    m_cascades( new QSpinBox )
{
    setWindowTitle(tr("Выбор параметров сглаживающего фильтра"));
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    setMinimumSize(400, 180);
    setToolTip("рабочие варианты: {2: 5..11}; {3: 5..11}; {4: 7..11}; {5: 7}");

    QGridLayout *layout = new QGridLayout(this);

    layout->addWidget(new QLabel(tr("порядок полинома :")), 0, 0);
    layout->addWidget(m_n_poly, 0, 1);
    m_n_poly->setRange( 1, 7 );
    m_n_poly->setValue( n_poly );
    m_n_poly->setAlignment( Qt::AlignCenter );

    layout->addWidget(new QLabel(tr("полуширина фильтра:")), 1, 0);
    layout->addWidget( m_halfwidth, 1, 1);
    m_halfwidth->setRange( 1, 9 );
    m_halfwidth->setValue( halfwidth );
    m_halfwidth->setAlignment( Qt::AlignCenter );

    layout->addWidget(new QLabel(tr("число каскадов:")), 2, 0);
    layout->addWidget( m_cascades, 2, 1);
    m_cascades->setRange( 1, 5 );
    m_cascades->setValue( cascades );
    m_cascades->setAlignment( Qt::AlignCenter );

    QDialogButtonBox *buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);
    connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
    layout->addWidget(buttonBox, 3, 0, 1, 2);
}

SFilterParams* DefineFilterParamDialog::getFilterParams() const
{
    SFilterParams* FilterParams = new SFilterParams;
    FilterParams->n_poly = m_n_poly->value();
    DefineFilterParamDialog::n_poly = FilterParams->n_poly;
    FilterParams->halfwidth = m_halfwidth->value();
    DefineFilterParamDialog::halfwidth = FilterParams->halfwidth;
    FilterParams->cascades = m_cascades->value();
    DefineFilterParamDialog::cascades = FilterParams->cascades;


    return FilterParams;
};
