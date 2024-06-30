#!/usr/bin/env python
# -*- coding:utf-8 -*-

# 本文件只允许依赖math库
import math


def draw_line(p_list, algorithm):
    """绘制线段

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'，此处的'Naive'仅作为示例，测试时不会出现
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    result = []
    if algorithm == 'Naive':
        if x0 == x1:
            for y in range(y0, y1 + 1):
                result.append((x0, y))
        else:
            if x0 > x1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            k = (y1 - y0) / (x1 - x0)
            for x in range(x0, x1 + 1):
                result.append((x, int(y0 + k * (x - x0))))
    elif algorithm == 'DDA':
        dx = x1 - x0
        dy = y1 - y0
        steps = max(abs(dx), abs(dy))
        if steps == 0:
            return result
        x_inc = dx / steps
        y_inc = dy / steps
        x, y = x0, y0
        for i in range(steps + 1):
            result.append([int(x), int(y)])
            x += x_inc
            y += y_inc
    elif algorithm == 'Bresenham':
        """NOTE: Bresenham General Line Drawing Algorithm"""
        dx, dy = abs(x1 - x0), abs(y1 - y0)
        sx = 1 if x0 < x1 else -1
        sy = 1 if y0 < y1 else -1
        err = dx - dy

        while True:
            result.append((x0, y0))
            if x0 == x1 and y0 == y1:
                break
            e2 = 2 * err
            if e2 > -dy:
                err -= dy
                x0 += sx
            if e2 < dx:
                err += dx
                y0 += sy
    return result


def draw_polygon(p_list, algorithm):
    """绘制多边形

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 多边形的顶点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    for i in range(len(p_list)):
        line = draw_line([p_list[i - 1], p_list[i]], algorithm)
        result += line
    return result


def draw_ellipse(p_list):
    """绘制椭圆（采用中点圆生成算法）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 椭圆的矩形包围框左上角和右下角顶点坐标
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    # NOTE: Midpoint Ellipse Drawing Algorithm
    x0, y0, x1, y1 = p_list[0][0], p_list[0][1], p_list[1][0], p_list[1][1]
    if x0 > x1:
        x0, x1 = x1, x0
    if y0 > y1:
        y0, y1 = y1, y0
    rx = (x1 - x0) // 2
    ry = (y1 - y0) // 2
    xc = x0 + rx
    yc = y0 + ry
    points = []

    x, y = 0, ry
    p1 = ry**2 - rx**2 * ry + 0.25 * rx**2
    while 2 * ry**2 * x < 2 * rx**2 * y:
        points += [[xc + x, yc + y], [xc - x, yc + y], [xc + x, yc - y], [xc - x, yc - y]]
        if p1 < 0:
            x += 1
            p1 += 2 * ry**2 * x + ry**2
        else:
            x += 1
            y -= 1
            p1 += 2 * ry**2 * x - 2 * rx**2 * y + ry**2

    p2 = ry**2 * (x + 0.5)**2 + rx**2 * (y - 1)**2 - rx**2 * ry**2
    while y >= 0:
        points += [[xc + x, yc + y], [xc - x, yc + y], [xc + x, yc - y], [xc - x, yc - y]]
        if p2 > 0:
            y -= 1
            p2 -= 2 * rx**2 * y + rx**2
        else:
            y -= 1
            x += 1
            p2 += 2 * ry**2 * x - 2 * rx**2 * y + rx**2

    return points

def bezier_curve(points, t):
    """计算给定时间t时贝塞尔曲线上的点"""
    n = len(points) - 1
    x = sum([points[i][0] * bernstein_poly(i, n, t) for i in range(n + 1)])
    y = sum([points[i][1] * bernstein_poly(i, n, t) for i in range(n + 1)])
    return [int(x), int(y)]

def bernstein_poly(i, n, t):
    """计算贝恩斯坦多项式的值"""
    return comb(n, i) * (t ** i) * ((1 - t) ** (n - i))


# 辅助函数，用于计算组合数，即从n个不同元素中取出k个元素的组合数
def comb(n, k):
    return factorial(n) / (factorial(k) * factorial(n - k))

# 阶乘函数
def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)

import numpy as np

def b_spline_basis(i, k, t, u):
    """计算B样条的基函数，k为基函数的度，i为控制点索引，t为参数，u为节点向量"""
    if k == 0:
        return 1.0 if u[i] <= t < u[i+1] else 0.0
    if u[i+k] == u[i]:
        c1 = 0.0
    else:
        c1 = (t - u[i]) / (u[i+k] - u[i]) * b_spline_basis(i, k-1, t, u)
    if u[i+k+1] == u[i+1]:
        c2 = 0.0
    else:
        c2 = (u[i+k+1] - t) / (u[i+k+1] - u[i+1]) * b_spline_basis(i+1, k-1, t, u)
    return c1 + c2

def draw_b_spline(points, u):
    """绘制B样条曲线"""
    n = len(points) - 1
    k = 3  # B样条曲线的度
    curve_points = []

    for t in np.linspace(u[k], u[n+1], num=1000):
        x, y = 0.0, 0.0
        for i in range(n+1):
            coeff = b_spline_basis(i, k, t, u)
            x += coeff * points[i][0]
            y += coeff * points[i][1]
        curve_points.append([int(x), int(y)])
    return curve_points


def draw_curve(p_list, algorithm):
    """绘制曲线

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 曲线的控制点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'Bezier'和'B-spline'（三次均匀B样条曲线，曲线不必经过首末控制点）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    curve_points = []
    if algorithm == 'Bezier':
        for t in [i / 1000.0 for i in range(1001)]:
            curve_points.append(bezier_curve(p_list, t))
    elif algorithm == 'B-spline':
        # 假设节点向量为均匀的情况
        n = len(p_list) + 3  # 控制点数加上度数（3）的和
        u = np.linspace(0, 1, n+1)  # 均匀分布的节点向量
        return draw_b_spline(p_list, u)
    else:
        raise ValueError("Unsupported algorithm")
    return curve_points


def translate(p_list, dx, dy):
    """平移变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param dx: (int) 水平方向平移量
    :param dy: (int) 垂直方向平移量
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    for i in range(len(p_list)):
        p_list[i][0] += dx
        p_list[i][1] += dy
    return p_list


def rotate(p_list, x, y, r):
    """旋转变换（除椭圆外）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 旋转中心x坐标
    :param y: (int) 旋转中心y坐标
    :param r: (int) 顺时针旋转角度（°）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    for i in range(len(p_list)):
        x0, y0 = p_list[i]
        x0 -= x
        y0 -= y
        x1 = x0 * math.cos(math.radians(r)) - y0 * math.sin(math.radians(r))
        y1 = x0 * math.sin(math.radians(r)) + y0 * math.cos(math.radians(r))
        p_list[i] = [int(x1 + x), int(y1 + y)]
    return p_list


def scale(p_list, x, y, s):
    """缩放变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 缩放中心x坐标
    :param y: (int) 缩放中心y坐标
    :param s: (float) 缩放倍数
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    for i in range(len(p_list)):
        x0, y0 = p_list[i]
        x0 -= x
        y0 -= y
        x1 = x0 * s
        y1 = y0 * s
        p_list[i] = [int(x1 + x), int(y1 + y)]
    return p_list


def compute_out_code(x, y, x_min, y_min, x_max, y_max):
    """计算点的OutCode."""
    code = 0
    if x < x_min:  # 点在裁剪窗口左侧
        code |= 1
    elif x > x_max:  # 点在裁剪窗口右侧
        code |= 2
    if y < y_min:  # 点在裁剪窗口下方
        code |= 4
    elif y > y_max:  # 点在裁剪窗口上方
        code |= 8
    return code

def cohen_sutherland_clip(p_list, x_min, y_min, x_max, y_max):
    """使用Cohen-Sutherland算法进行线段裁剪."""
    x0, y0, x1, y1 = p_list[0][0], p_list[0][1], p_list[1][0], p_list[1][1]
    out_code0 = compute_out_code(x0, y0, x_min, y_min, x_max, y_max)
    out_code1 = compute_out_code(x1, y1, x_min, y_min, x_max, y_max)
    accept = False

    while True:
        if not (out_code0 | out_code1):
            accept = True
            break
        elif out_code0 & out_code1:
            break
        else:
            if out_code0:
                code_out = out_code0
            else:
                code_out = out_code1

            if code_out & 8:
                x = x0 + (x1 - x0) * (y_max - y0) / (y1 - y0)
                y = y_max
            elif code_out & 4:
                x = x0 + (x1 - x0) * (y_min - y0) / (y1 - y0)
                y = y_min
            elif code_out & 2:
                y = y0 + (y1 - y0) * (x_max - x0) / (x1 - x0)
                x = x_max
            elif code_out & 1:
                y = y0 + (y1 - y0) * (x_min - x0) / (x1 - x0)
                x = x_min

            if code_out == out_code0:
                x0, y0 = x, y
                out_code0 = compute_out_code(x0, y0, x_min, y_min, x_max, y_max)
            else:
                x1, y1 = x, y
                out_code1 = compute_out_code(x1, y1, x_min, y_min, x_max, y_max)

    if accept:
        return [[int(x0), int(y0)], [int(x1), int(y1)]]
    else:
        return []  # 裁剪窗口外的线段不返回任何点



def clip(p_list, x_min, y_min, x_max, y_max, algorithm):
    """线段裁剪

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param x_min: 裁剪窗口左上角x坐标
    :param y_min: 裁剪窗口左上角y坐标
    :param x_max: 裁剪窗口右下角x坐标
    :param y_max: 裁剪窗口右下角y坐标
    :param algorithm: (string) 使用的裁剪算法，包括'Cohen-Sutherland'和'Liang-Barsky'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1]]) 裁剪后线段的起点和终点坐标
    """
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    if algorithm == 'Cohen-Sutherland':
        # NOTE: Cohen-Sutherland Algorithm
        return cohen_sutherland_clip(p_list, x_min, y_min, x_max, y_max)
    elif algorithm == 'Liang-Barsky':
        # TODO: Implement Liang-Barsky Algorithm
        return Liang_Barsky_clip(p_list, x_min, y_min, x_max, y_max)
    return [[x0, y0], [x1, y1]]


def Liang_Barsky_clip(p_list, x_min, y_min, x_max, y_max):
    """Liang-Barsky 线段裁剪

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param x_min: 裁剪窗口左上角x坐标
    :param y_min: 裁剪窗口左上角y坐标
    :param x_max: 裁剪窗口右下角x坐标
    :param y_max: 裁剪窗口右下角y坐标
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1]]) 裁剪后线段的起点和终点坐标
    """
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]

    dx = x1 - x0
    dy = y1 - y0

    p = [-dx, dx, -dy, dy]
    q = [x0 - x_min, x_max - x0, y0 - y_min, y_max - y0]

    u1 = 0.0
    u2 = 1.0

    for i in range(4):
        if p[i] == 0:
            if q[i] < 0:
                return []  # Line is parallel to and outside the clipping window
        else:
            t = q[i] / p[i]
            if p[i] < 0:
                if t > u2:
                    return []  # Line is outside the clipping window
                elif t > u1:
                    u1 = t
            else:
                if t < u1:
                    return []  # Line is outside the clipping window
                elif t < u2:
                    u2 = t

    clipped_x0 = int(x0 + u1 * dx)
    clipped_y0 = int(y0 + u1 * dy)
    clipped_x1 = int(x0 + u2 * dx)
    clipped_y1 = int(y0 + u2 * dy)

    return [[clipped_x0, clipped_y0], [clipped_x1, clipped_y1]]