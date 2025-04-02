def print_dict(topDict):
    for key1, value1 in topDict.items():
        print(f"--> KEY 1\t {key1}")
        if type(value1) == dict:
            for key2, value2 in value1.items():
                print(f"{'-->'*2} KEY 2\t {key2}")
                if type(value2) == dict:
                    for key3, value3 in value2.items():
                        print(f"{'-->'*3} KEY 3\t {key3}")
                        print(f"{'-->'*3} VALUE 3\t {value3}")
                elif type(value2) == list:
                    for value3 in value2:
                        print(f"{'-->'*3} VALUE 2\t {value3}")
                else:
                    print(f"{'-->'*2} VALUE 2\t{value2}")

        elif type(value1) == list:
            for value2 in value1:
                print(f"{'-->'*2} VALUE 1\t{value2}")
        else:
            print(f"{'-->'*1} VALUE 1\t{value1}")


