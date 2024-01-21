import math
import time
import decimal
from fractions import Fraction
# decimal.Decimal

class HammingSeriesSum:
    def __init__(self) -> None:
        pass

    def HammingSeriesSumRecursiveSolver(self, x):
        x, sum, k = decimal.Decimal(x), decimal.Decimal("0.0"), decimal.Decimal("1.0")
        while k <= decimal.Decimal("60000.0"):
            sum += decimal.Decimal("1.0")/(k*(k+decimal.Decimal("1.0"))*(k+decimal.Decimal("2.0"))*(k+x))
            k += decimal.Decimal("1.0")
        sum = decimal.Decimal("1.0") + (decimal.Decimal("1.0")-x) * (decimal.Decimal("2.0")-x)* sum + (decimal.Decimal("1.0")-x)/decimal.Decimal("4.0")
        print(sum)
        return None


    def compare(self):
        x, sum, k = 0.0, 0.0, 1.0
        while k <= 1494.0:
            sum += 1.0/(k*(k+1.0)*(k+2.0)*(k+x))
            k += 1.0
        sum = 1.0 + (1.0-x) * ((2.0-x)* sum + 0.25)
        print(sum,k)

def HammingSeriesSumRecursiveSolver(self, x):
    sum, k = 0.0, 1.0
    while k <= 60000.0:
        sum += 1.0/(k*(k+1.0)*(k+2.0)*(k+x))
        k += 1.0
    sum = 1.0 + (1.0-x) * (2.0-x)* sum + (1.0-x)/4.0
    print(sum)
    return None


if __name__ == '__main__':

    hss = HammingSeriesSum()

    start_time = time.time()
    x = 0.0
    while x<=300.0:
        hss.fierce(x)
        x += 0.1
    end_time = time.time()
    print("time is {} ms".format((end_time-start_time)*1000))
    
    start_time = time.time()
    hss.compare()
    end_time = time.time()
    print("compare time is {} ms".format((end_time-start_time)*1000))