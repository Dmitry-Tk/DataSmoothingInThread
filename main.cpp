#include <QObject>
#include "ccontrolwidget.h"

/**
 * Начало: 12.04.2018
 *
 *      Полиномиальное сглаживание непрерывного потока данных с помощью цифрового фильтра.
 *      Коэфф-ты цифрового фильтра вычисляются автоматически при выборе порядка сглаживающего
 *      полинома и количества используемых соседних точек.
 *      Цифровой фильтр вычисляется исходя из критерия минимизации стандартного отклонения для
 *      сглаживающего полинома
 *
 *  1. +a) Организовать конвейерные вычисления в нити
 *     +b) Измерить скорость обработки данных в нити и оптимизировать её, если возможно
 *     -c) Сравнить со скоростью обработки потоков данных в старой программе:
 *               DemoPlotterThreadWorkerClass
 *  +2. Разработать генератор произвольного фильтра на min(STD) для класса фильтра
 *   3. Разработать класс интерполятора
 *  +4. Прекращение выполнения нити после закрытия программы
 *  +5. Добавить вычисление StD
**/

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    QMainWindow MainWindow;
    QDockWidget *DockCtrlWgt = new QDockWidget;
    CControlWidget *ControlWidget = new CControlWidget;
    ///
    /// виджет управления
    DockCtrlWgt->setWidget(ControlWidget);
    DockCtrlWgt->setFeatures( QDockWidget::DockWidgetMovable );
//    DockCtrlWgt->setFeatures( QDockWidget::DockWidgetFloatable );
    MainWindow.addDockWidget(Qt::TopDockWidgetArea, DockCtrlWgt);
    /// инсталлируем виджет плоттера
    MainWindow.setCentralWidget( ControlWidget->ChartView );
    /// главное окно
    MainWindow.setWindowTitle("Полиномиальное сглаживание (фильтрация) данных минимизацией STD");
    MainWindow.setGeometry(100, 100, 1100, 700);
    MainWindow.show();

    return app.exec();
}
