#!/usr/bin/env python3
#coding=utf-8

from PySpice.Unit.Units import *

E12R=[1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1, 10, 11, 12, 13, 15, 16, 18, 20, 22, 24, 27, 30, 33, 36, 39, 43, 47, 51, 56, 62, 68, 75, 82, 91, 100, 110, 120, 130, 150, 160, 180, 200, 220, 240, 270, 300, 330, 360, 390, 430, 470, 510, 560, 620, 680, 750, 820, 910, 1.0, kilo(1.1), kilo(1.2), kilo(1.3), kilo(1.5), kilo(1.6), kilo(1.8), kilo(2.0), kilo(2.2), kilo(2.4), kilo(2.7), kilo(3.0), kilo(3.3), kilo(3.6), kilo(3.9), kilo(4.3), kilo(4.7), kilo(5.1), kilo(5.6), kilo(6.2), kilo(6.8), kilo(7.5), kilo(8.2), kilo(9.1), kilo(10), kilo(11), kilo(12), kilo(13), kilo(15), kilo(16), kilo(18), kilo(20), kilo(22), kilo(24), kilo(27), kilo(30), kilo(33), kilo(36), kilo(39), kilo(43), kilo(47), kilo(51), kilo(56), kilo(62), kilo(68), kilo(75), kilo(82), kilo(91), kilo(100), kilo(110), kilo(120), kilo(130), kilo(150), kilo(160), kilo(180), kilo(200), kilo(220), kilo(240), kilo(270), kilo(300), kilo(330), kilo(360), kilo(390), kilo(430), kilo(470), kilo(510), kilo(5-60), kilo(620), kilo(680), kilo(750), kilo(820), kilo(910), mega(1.0) ,mega(1.1) ,mega(1.2) ,mega(1.3) ,mega(1.5) ,mega(1.6) ,mega(1.8) ,mega(2.0) ,mega(2.2) ,mega(2.4) ,mega(2.7) ,mega(3.0) ,mega(3.3) ,mega(3.6) ,mega(3.9) ,mega(4.3) ,mega(4.7) ,mega(5.1) ,mega(5.6) ,mega(6.2) ,mega(6.8) ,mega(7.5) ,mega(8.2) ,mega(9.1)]

E12C=[pico(10), pico(12), pico(15), pico(18), pico(22), pico(27), pico(33), pico(39), pico(47), pico(56), pico(68), pico(82), pico(100), pico(120), pico(150), pico(180), pico(220), pico(270), pico(330), pico(390), pico(470), pico(560), pico(680), pico(820), pico(1000), pico(1200), pico(1500), pico(1800), pico(2200), pico(2700), pico(3300), pico(3900), pico(4700), pico(5600), pico(6800), pico(8200), micro(.010), micro(.012), micro(.015), micro(.018), micro(.022), micro(.027), micro(.033), micro(.039), micro(.047), micro(.056), micro(.068), micro(.082), micro(.10), micro(.12), micro(.15), micro(.18), micro(.22), micro(.27), micro(.33), micro(.39), micro(.47), micro(.56), micro(.68), micro(.82), micro(1.0), micro(1.2), micro(1.5), micro(1.8), micro(2.2), micro(2.7), micro(3.3), micro(3.9), micro(4.7), micro(5.6), micro(6.8), micro(8.2), micro(10), micro(22), micro(33), micro(47)]

E12I=[micro(.010), micro(.012), micro(.015), micro(.018), micro(.022), micro(.027), micro(.033), micro(.039), micro(.047), micro(.056), micro(.068), micro(.082), micro(.10), micro(.12), micro(.15), micro(.18), micro(.22), micro(.27), micro(.33), micro(.39), micro(.47), micro(.56), micro(.68), micro(.82), micro(1.0), micro(1.2), micro(1.5), micro(1.8), micro(2.2), micro(2.7), micro(3.3), micro(3.9), micro(4.7), micro(5.6), micro(6.8), micro(8.2), micro(10), micro(22), micro(33), micro(47), micro(100), micro(120), micro(150), micro(180), micro(220), micro(270), micro(330), micro(390), micro(470), micro(560), micro(680), micro(820), micro(100), micro(220), micro(330), micro(470)]
