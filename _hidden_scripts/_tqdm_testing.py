from tqdm import tqdm
from time import sleep
tqdmBarOptions = {
    "desc": f"\033[32mRunning Single-Points \033[0m",
    "ascii": "-🗲→",  
    "colour": "yellow",
    "unit":  "scan",
    "dynamic_ncols": True,
    "leave": True,

}


totalTickProgress =  100.00
progress_bar = tqdm(total=totalTickProgress, **tqdmBarOptions)

# for _ in range(100):
#     progress_bar.update(1)
#     sleep(0.05)


class ColoredTqdm(tqdm):
    def format_meter(self, n, total, elapsed, **kwargs):
        greenText = "\033[32m"
        redText = "\033[31m"
        orangeText = "\033[38;5;172m"
        yellowText = "\033[33m"
        resetTextColor = "\033[0m"
        purpleText = "\033[35m"

        # Calculate bar width dynamically, accounting for desc and other text
        desc_len = len(self.desc or '') + 2  # +2 for padding
        bar_width = (self.ncols or 50) - desc_len - 15  # 15 for "n/total [time]"
        filled = int(bar_width * n / total)
        unfilled = bar_width - filled

        # Handle edge cases
        if filled <= 0:
            bar = f"{purpleText}{'-' * bar_width}{resetTextColor}"
        elif filled >= bar_width:
            bar = f"{greenText}{'→' * (bar_width - 1)}{yellowText}{'🗲'}{resetTextColor}"
        else:
            bar = f"{greenText}{'→' * (filled - 1)}{yellowText}{'🗲'}{purpleText}{'-' * unfilled}{resetTextColor}"

        return f"{self.desc or ''} {bar} {n}/{total} [{elapsed:.1f}s]"

# Usage
for i in ColoredTqdm(range(100), desc="\033[33mScan\033[0m", unit="step", dynamic_ncols=True):
    sleep(0.05)