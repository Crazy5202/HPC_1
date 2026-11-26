import matplotlib.pyplot as plt

# Данные
x = [1, 2, 4]

y_100 = [1.0, 0.000956/0.000706, 0.000956/0.000657]
y_1000 = [1.0, 0.105401/0.061207, 0.105401/0.038297]
y_10000 = [1.0, 11.458012/6.230278, 11.458012/3.640418]

# Создаем график
plt.plot(x, x, marker='o', linestyle='--', color='black', label = 'Идеальное ускорение')
plt.plot(x, y_100, marker='o', linestyle='-', color='r', label = 'Ускорение при m=100')
plt.plot(x, y_1000, marker='o', linestyle='-', color='b', label = 'Ускорение при m=1 000')
plt.plot(x, y_10000, marker='o', linestyle='-', color='g', label = 'Ускорение при m=10 000')

# Настройки осей и заголовок
plt.xlabel('Число задействованных потоков')
plt.ylabel('Ускорение')

plt.xticks(x)
plt.yticks(x+y_100[1:]+y_1000[1:]+y_10000[1:])

plt.legend()

plt.title('Ускорение в зависимости от числа потоков')

plt.grid(True, which="both", ls="-", alpha=0.2)

# Отображаем график
plt.show()