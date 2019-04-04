import sys
sys.path.append("../../")
from driver import GaussDriver, ExitSignal, TimeOutSignal

driver = GaussDriver(31415, "127.0.0.1", "template.gjf", ["H", "H"])
while True:
    try:
        driver.parse()
    except ExitSignal as e:
        driver = GaussDriver(31415, "127.0.0.1", "template.gjf", ["H", "H"])
    except TimeOutSignal as e:
        print("Time out. Check whether the server is closed.")
        exit()