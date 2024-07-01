#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import cg_algorithms as alg
import math
from typing import Optional
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    qApp,
    QGraphicsScene,
    QGraphicsView,
    QGraphicsItem,
    QListWidget,
    QHBoxLayout,
    QWidget,
    QGraphicsRectItem,
    QColorDialog,
    QInputDialog,
    QStyleOptionGraphicsItem)
from PyQt5.QtGui import QPainter, QMouseEvent, QColor
from PyQt5.QtCore import QRectF, Qt


class MyCanvas(QGraphicsView):
    """
    画布窗体类，继承自QGraphicsView，采用QGraphicsView、QGraphicsScene、QGraphicsItem的绘图框架
    """
    def __init__(self, *args):
        super().__init__(*args)
        self.main_window = None
        self.list_widget = None
        self.item_dict = {}
        self.selected_id = ''

        self.status = ''
        self.temp_algorithm = ''
        self.temp_id = ''
        self.temp_item = None

        self.polygon = 0
        self.curve = 0

    def start_draw_line(self, algorithm, item_id):
        self.status = 'line'
        self.temp_algorithm = algorithm
        self.temp_id = item_id

    def start_draw_polygon(self, algorithm, item_id):
        self.status = 'polygon'
        self.temp_algorithm = algorithm
        self.temp_id = item_id
    
    def start_draw_ellipse(self, algorithm, item_id):
        self.status = 'ellipse'
        self.temp_algorithm = algorithm
        self.temp_id = item_id 

    def start_draw_curve(self, algorithm, item_id):
        self.status = 'curve'
        self.temp_algorithm = algorithm
        self.temp_id = item_id
        
    def finish_draw(self):
        self.temp_id = self.main_window.get_id()

    def start_translate(self):
        self.status = 'translate'

    def start_rotate(self):
        self.status = 'rotate'

    def start_scale(self):
        self.status = 'scale'

    def start_clip(self, algorithm):
        self.status = 'clip'
        self.temp_algorithm = algorithm
        self.clipRect = None

    def clear_selection(self):
        if self.selected_id != '':
            self.item_dict[self.selected_id].selected = False
            self.selected_id = ''

    def selection_changed(self, selected):
        if selected != '':
            self.main_window.statusBar().showMessage('图元选择： %s' % selected)
            if self.selected_id != '':
                self.item_dict[self.selected_id].selected = False
                self.item_dict[self.selected_id].update()
            self.selected_id = selected
            self.item_dict[selected].selected = True
            self.item_dict[selected].update()
            self.status = ''
            self.updateScene([self.sceneRect()])

    def mousePressEvent(self, event: QMouseEvent) -> None:
        pos = self.mapToScene(event.localPos().toPoint())
        x = int(pos.x())
        y = int(pos.y())
        if self.status == 'line' or self.status == 'ellipse':
            self.temp_item = MyItem(self.temp_id, self.status, [[x, y], [x, y]], self.temp_algorithm)
            self.scene().addItem(self.temp_item)
        elif self.status == 'polygon':
            if event.buttons() == Qt.RightButton and self.polygon == 1:
                self.list_widget.addItem(self.temp_id)
                self.finish_draw()
                self.polygon = 0
            elif event.buttons() == Qt.LeftButton:
                if self.polygon == 0:
                    self.polygon = 1
                    self.temp_item = MyItem(self.temp_id, self.status, [[x, y], [x, y]], self.temp_algorithm) 
                    self.scene().addItem(self.temp_item)
                else:
                    self.temp_item.p_list.append([x, y])
                
        elif self.status == 'curve':
            if event.buttons() == Qt.RightButton and self.curve == 1:
                self.list_widget.addItem(self.temp_id)
                self.finish_draw()
                self.curve = 0
            elif event.buttons() == Qt.LeftButton:
                if self.curve == 0:
                    self.curve = 1
                    self.temp_item = MyItem(self.temp_id, self.status, [[x, y], [x, y]], self.temp_algorithm) 
                    self.scene().addItem(self.temp_item)
                else:
                    self.temp_item.p_list.append([x, y])
            
        elif self.status in ['translate', 'rotate', 'scale']:
            if self.selected_id != '':
                self.temp_item = self.item_dict[self.selected_id]
                self.x_old, self.y_old = x, y
                self.p_list_old = self.temp_item.p_list
        elif self.status == 'clip':
            self.temp_item = None
            if self.selected_id != '' and self.item_dict[self.selected_id] != None:
                if self.item_dict[self.selected_id].item_type == 'line': 
                    self.temp_item = self.item_dict[self.selected_id]
                    self.x_old, self.y_old = x, y
                    self.p_list_old = self.temp_item.p_list
        self.updateScene([self.sceneRect()])
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        pos = self.mapToScene(event.localPos().toPoint())
        x = int(pos.x())
        y = int(pos.y())
        if self.status == 'line' or self.status == 'ellipse':
            self.temp_item.p_list[1] = [x, y]
        elif self.status == 'polygon':
            self.temp_item.p_list[-1] = [x, y]
        elif self.status == 'curve':
            self.temp_item.p_list[-1] = [x, y]
        elif self.status == 'translate':
            if self.selected_id != '':
                dx, dy = x - self.x_old, y - self.y_old
                self.temp_item.p_list = alg.translate(self.p_list_old, dx, dy)
        elif self.status == 'rotate':
            if self.selected_id != '':
                dx, dy = x - self.x_old, y - self.y_old
                rad = math.atan2(dy, dx)
                angle = math.degrees(rad)
                self.temp_item.p_list = alg.rotate(self.p_list_old, self.x_old, self.y_old, angle)
        elif self.status == 'scale':
            if self.selected_id != '':
                s = y / max(1., self.y_old)
                self.temp_item.p_list = alg.scale(self.p_list_old, self.x_old, self.y_old, s)
        elif self.status == 'clip':
            if self.temp_item != None and self.temp_item.item_type == 'line':
                x_min = min(x, self.x_old)
                x_max = max(x, self.x_old)
                y_min = min(y, self.y_old)
                y_max = max(y, self.y_old)
                w = x_max - x_min
                h = y_max - y_min
                if self.clipRect is None:
                    self.clipRect = QGraphicsRectItem(x_min - 1, y_min - 1, w+2, h+2)
                    self.scene().addItem(self.clipRect)
                    self.clipRect.setPen(QColor(0, 233, 255))
                else:
                    self.clipRect.setRect(x_min - 1, y_min - 1, w + 2, h + 2)

        self.updateScene([self.sceneRect()])
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        if self.status == 'line' or self.status == 'ellipse':
            self.item_dict[self.temp_id] = self.temp_item
            self.list_widget.addItem(self.temp_id)
            self.finish_draw()
        elif self.status == 'polygon':
            self.item_dict[self.temp_id] = self.temp_item
        elif self.status == 'curve':
            self.item_dict[self.temp_id] = self.temp_item
        elif self.status == 'clip':
            pos = self.mapToScene(event.localPos().toPoint())
            x = int(pos.x())
            y = int(pos.y())
            if self.temp_item != None and self.temp_item.item_type=='line':
                x_min = min(self.x_old, x)
                x_max = max(self.x_old, x)
                y_min = min(self.y_old, y)
                y_max = max(self.y_old, y)
                temp_p_list = alg.clip(self.p_list_old, x_min, y_min, x_max, y_max, self.temp_algorithm)
                if len(temp_p_list) == 0:
                    number = self.list_widget.findItems(self.selected_id, Qt.MatchContains)
                    row = self.list_widget.row(number[0])
                    temp_id = self.selected_id
                    self.clear_selection()
                    self.list_widget.clearSelection()
                    self.scene().removeItem(self.temp_item)
                    self.temp_item = None
                    del self.item_dict[temp_id]
                    self.list_widget.takeItem(row)
                    # self.temp_item.p_list = temp_p_list
                else:
                    self.temp_item.p_list = temp_p_list
                if self.clipRect is not None:
                    self.scene().removeItem(self.clipRect)
                    self.clipRect = None
                self.updateScene([self.sceneRect()])
        super().mouseReleaseEvent(event)


class MyItem(QGraphicsItem):
    """
    自定义图元类，继承自QGraphicsItem
    """
    def __init__(self, item_id: str, item_type: str, p_list: list, algorithm: str = '', parent: QGraphicsItem = None):
        """

        :param item_id: 图元ID
        :param item_type: 图元类型，'line'、'polygon'、'ellipse'、'curve'等
        :param p_list: 图元参数
        :param algorithm: 绘制算法，'DDA'、'Bresenham'、'Bezier'、'B-spline'等
        :param parent:
        """
        super().__init__(parent)
        self.id = item_id           # 图元ID
        self.item_type = item_type  # 图元类型，'line'、'polygon'、'ellipse'、'curve'等
        self.p_list = p_list        # 图元参数
        self.algorithm = algorithm  # 绘制算法，'DDA'、'Bresenham'、'Bezier'、'B-spline'等
        self.selected = False
        self.color = mw.canvas_widget.color

    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[QWidget] = ...) -> None:
        painter.setPen(self.color)
        if self.item_type == 'line':
            item_pixels = alg.draw_line(self.p_list, self.algorithm)
            for p in item_pixels:
                painter.drawPoint(*p)
            if self.selected:
                painter.setPen(QColor(255, 0, 0))
                painter.drawRect(self.boundingRect())
        elif self.item_type == 'polygon':
            item_pixels = alg.draw_polygon(self.p_list, self.algorithm)
            # print("item_pixels:\n", item_pixels)
            for p in item_pixels:
                painter.drawPoint(*p)
            if self.selected:
                painter.setPen(QColor(0, 0, 255))
                painter.drawRect(self.boundingRect()) 
        elif self.item_type == 'ellipse':
            item_pixels = alg.draw_ellipse(self.p_list)
            for p in item_pixels:
                painter.drawPoint(*p)
            if self.selected:
                painter.setPen(QColor(0, 255, 0))
                painter.drawRect(self.boundingRect())
        elif self.item_type == 'curve':
            if self.algorithm == 'Bezier' or len(self.p_list) > 3:
                item_pixels = alg.draw_curve(self.p_list, self.algorithm)
                for p in item_pixels:
                    painter.drawPoint(*p)
                if self.selected:
                    painter.setPen(QColor(255, 0, 255))
                    painter.drawRect(self.boundingRect())

    def boundingRect(self) -> QRectF:
        if self.item_type == 'line' or self.item_type == 'ellipse':
            x0, y0 = self.p_list[0]
            x1, y1 = self.p_list[1]
            x = min(x0, x1)
            y = min(y0, y1)
            w = max(x0, x1) - x
            h = max(y0, y1) - y
            return QRectF(x - 1, y - 1, w + 2, h + 2)
        elif self.item_type == 'polygon' or self.item_type == 'curve':
            x_min, y_min = self.p_list[0]
            x_max, y_max = self.p_list[0]
            for p in self.p_list:
                x_min = min(x_min, p[0])
                x_max = max(x_max, p[0])
                y_min = min(y_min, p[1])
                y_max = max(y_max, p[1])
            w = x_max - x_min
            h = y_max - y_min
            return QRectF(max(x_min - 1, 0), max(y_min - 1,0), w + 2, h + 2)
        else:
            return QRectF(0, 0, 600, 600) 


class MainWindow(QMainWindow):
    """
    主窗口类
    """
    def __init__(self):
        super().__init__()
        self.item_cnt = 0

        # 使用QListWidget来记录已有的图元，并用于选择图元。注：这是图元选择的简单实现方法，更好的实现是在画布中直接用鼠标选择图元
        self.list_widget = QListWidget(self)
        self.list_widget.setMinimumWidth(200)

        # 使用QGraphicsView作为画布
        self.scene = QGraphicsScene(self)
        self.scene.setSceneRect(0, 0, 600, 600)
        self.canvas_widget = MyCanvas(self.scene, self)
        self.canvas_widget.setFixedSize(600, 600)
        self.canvas_widget.main_window = self
        self.canvas_widget.list_widget = self.list_widget
        self.canvas_widget.color = QColor(0, 0, 0)

        # 设置菜单栏
        menubar = self.menuBar()
        file_menu = menubar.addMenu('文件')
        set_pen_act = file_menu.addAction('设置画笔')
        reset_canvas_act = file_menu.addAction('重置画布')
        exit_act = file_menu.addAction('退出')
        draw_menu = menubar.addMenu('绘制')
        line_menu = draw_menu.addMenu('线段')
        line_naive_act = line_menu.addAction('Naive')
        line_dda_act = line_menu.addAction('DDA')
        line_bresenham_act = line_menu.addAction('Bresenham')
        polygon_menu = draw_menu.addMenu('多边形')
        polygon_dda_act = polygon_menu.addAction('DDA')
        polygon_bresenham_act = polygon_menu.addAction('Bresenham')
        ellipse_act = draw_menu.addAction('椭圆')
        curve_menu = draw_menu.addMenu('曲线')
        curve_bezier_act = curve_menu.addAction('Bezier')
        curve_b_spline_act = curve_menu.addAction('B-spline')
        edit_menu = menubar.addMenu('编辑')
        translate_act = edit_menu.addAction('平移')
        rotate_act = edit_menu.addAction('旋转')
        scale_act = edit_menu.addAction('缩放')
        clip_menu = edit_menu.addMenu('裁剪')
        clip_cohen_sutherland_act = clip_menu.addAction('Cohen-Sutherland')
        clip_liang_barsky_act = clip_menu.addAction('Liang-Barsky')

        # 连接信号和槽函数
        exit_act.triggered.connect(qApp.quit)
        set_pen_act.triggered.connect(self.set_pen_action)
        line_naive_act.triggered.connect(self.line_naive_action)
        line_dda_act.triggered.connect(self.line_DDA_action)
        line_bresenham_act.triggered.connect(self.line_bresenham_action)
        polygon_dda_act.triggered.connect(self.polygon_dda_action)
        polygon_bresenham_act.triggered.connect(self.polygon_bresenham_action)
        ellipse_act.triggered.connect(self.ellipse_action)
        curve_bezier_act.triggered.connect(self.curve_bezier_act)
        curve_b_spline_act.triggered.connect(self.curve_spline_act)
        translate_act.triggered.connect(self.translate_action)
        rotate_act.triggered.connect(self.rotate_action)
        scale_act.triggered.connect(self.scale_action)
        clip_cohen_sutherland_act.triggered.connect(self.clip_cohen_sutherland_action)
        clip_liang_barsky_act.triggered.connect(self.clip_liang_barsky_actiton)
        reset_canvas_act.triggered.connect(self.reset_canvas_action)
        self.list_widget.currentTextChanged.connect(self.canvas_widget.selection_changed)

        # 设置主窗口的布局
        self.hbox_layout = QHBoxLayout()
        self.hbox_layout.addWidget(self.canvas_widget)
        self.hbox_layout.addWidget(self.list_widget, stretch=1)
        self.central_widget = QWidget()
        self.central_widget.setLayout(self.hbox_layout)
        self.setCentralWidget(self.central_widget)
        self.statusBar().showMessage('空闲')
        self.resize(600, 600)
        self.setWindowTitle('CG Demo')

    def get_id(self):
        _id = str(self.item_cnt)
        self.item_cnt += 1
        return _id

    def set_pen_action(self):
        self.statusBar().showMessage('设置颜色')
        color = QColorDialog.getColor()
        if color.isValid():
            self.canvas_widget.color = color
            self.statusBar().showMessage(f"画笔颜色设置为 {color.name()}")
    
    def reset_canvas_action(self):
        pass

    def line_naive_action(self):
        self.canvas_widget.start_draw_line('Naive', self.get_id())
        self.statusBar().showMessage('Naive算法绘制线段')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection()

    def line_DDA_action(self):
        self.canvas_widget.start_draw_line('DDA', self.get_id())
        self.statusBar().showMessage('DDA算法绘制线段')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection()
    
    def line_bresenham_action(self):
        self.canvas_widget.start_draw_line('Bresenham', self.get_id())
        self.statusBar().showMessage('Bresenham算法绘制线段')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection() 

    def polygon_dda_action(self):
        self.canvas_widget.start_draw_polygon('DDA', self.get_id())
        self.statusBar().showMessage('DDA 算法绘制多边形')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection() 
        
    def polygon_bresenham_action(self):
        self.canvas_widget.start_draw_polygon('Bresenham', self.get_id())
        self.statusBar().showMessage('Bresenham 算法绘制多边形')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection() 

    def ellipse_action(self):
        self.canvas_widget.start_draw_ellipse('Nothing', self.get_id())
        self.statusBar().showMessage('绘制椭圆')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection() 
    
    def curve_bezier_act(self):
        self.canvas_widget.start_draw_curve('Bezier', self.get_id())
        self.statusBar().showMessage('Bezier 算法绘制曲线')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection()  

    def curve_spline_act(self):
        self.canvas_widget.start_draw_curve('B-spline', self.get_id())
        self.statusBar().showMessage('B-spline 算法绘制曲线')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection()  
    
    def translate_action(self):
        self.canvas_widget.start_translate()
        self.statusBar().showMessage('Translating')
    
    def rotate_action(self):
        self.canvas_widget.start_rotate()
        self.statusBar().showMessage('Rotating')
    
    def scale_action(self):
        self.canvas_widget.start_scale()
        self.statusBar().showMessage('Scaling')
    
    def clip_cohen_sutherland_action(self):
        self.canvas_widget.start_clip('Cohen-Sutherland')
        self.statusBar().showMessage('Cohen-Sutherland 算法裁剪线段')

    def clip_liang_barsky_actiton(self):
        self.canvas_widget.start_clip('Liang-Barsky')
        self.statusBar().showMessage('Liang-Barsky 算法裁剪线段')
    
    

if __name__ == '__main__':
    app = QApplication(sys.argv)
    mw = MainWindow()
    mw.show()
    sys.exit(app.exec_())
