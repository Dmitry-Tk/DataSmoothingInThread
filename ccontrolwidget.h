#ifndef CCONTROLWIDGET_H
#define CCONTROLWIDGET_H
#include <QtWidgets>
#include <QtCharts/QChartView>
#include <QtCharts/QChart>
#include <QtCharts/QLineSeries>
#include <QtCharts/QScatterSeries>
#include <QPointF>
#include "ccalc.h"

using namespace QtCharts;

class CControlWidget : public QWidget
{
    Q_OBJECT

    QThread CalcThread;// ??? здесь ли

public:
    CControlWidget(QWidget *parent = 0);
    ~CControlWidget();
    /// пользовательский  интерфейс
    QPushButton* btnStartStop;
    QPushButton* btnSmoothData;
    QPushButton* btnResetChart;
    QPushButton* btnSetOriginalZoom;
    QPushButton* btnSetFilterParams;
    QVector<QRadioButton*> radioBtns;
    QSpinBox* spinboxNblocks;
    QSpinBox* spinboxDelay;
    QCheckBox* checkboxYScale;
    QSlider* sliderNoise;
    QChartView* ChartView;
    QChart* chart;
    // 0 - чистый; 1 - зашумленный; 2 - отфильтрованный сигналы
    QXYSeries *series[3];
    QVector<QPointF>* data[3];// вектора сигналов
    /// модель - сглаживающий непрерывный поток данных объект
    CCalc SmoothingObject;
    CBounds* Bounds;
    bool LogScale = false;

    void InitPlotter();
    void ConnectDataChannel();
    void CreateSeries();

public slots:
    void setSimulationMode();
    void setNoiseAmp(int);
    void setFilterParams();
    void StartStopCalc();
    void onUpdatePlotter();
    void onProcessFinished();
    void setYScale(bool);
    void ClearChart();
    void SetOriginalZoom();
    void SmoothData();
    void SmoothSeriesI( std::vector<QPointF>& Data, std::vector<QPointF>& SmoothedData);
    void SmoothSeriesII(const std::vector<QPointF>& Data, std::vector<QPointF>& SmoothedData);
    void SmoothSeriesIII(const std::vector<QPointF>& Data, std::vector<QPointF>& SmoothedData);

    void LoadData(const QString& FullFileNameString, std::vector<QPointF>& data);
    void SaveData(const QString& FullFileNameString, const std::vector<QPointF>& Data);


signals:
    void SimulationModeChanged(int);
};

struct SFilterParams
{
    int n_poly;/**< порядок сглаживающего полинома */
    int halfwidth;/**< полуширина фильтра в точках */
    int cascades;/**< число каскадов */
    QString m_name;/**< название фильтра */
};
///
/// \brief The DefineFilterParamDialog class
///
class DefineFilterParamDialog : public QDialog
{
public:
    explicit DefineFilterParamDialog(QWidget *parent = nullptr);
    SFilterParams* getFilterParams() const;

private:
    QSpinBox* m_n_poly;
    QSpinBox* m_halfwidth;
    QSpinBox* m_cascades;
    QString m_name;
    static int n_poly;/**< порядок сглаживающего полинома */
    static int halfwidth;/**< полуширина фильтра в точках */
    static int cascades;/**< число каскадов */
//    QString name;/**< название фильтра */
};

#endif // CCONTROLWIDGET_H
