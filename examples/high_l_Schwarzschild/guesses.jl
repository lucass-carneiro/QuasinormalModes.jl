# ------------------------------------------------------------------
# Reference QNMs for use as initial guesses
#
# Format: (l, n, real part, imaginary part)
# ------------------------------------------------------------------

const spin_0 = [
    (0, "0", "0.1104", "-0.1048"),
    (1, "0", "0.2911", "-0.0980"),
    (2, "0", "0.4832", "-0.0968"),
    (2, "1", "0.4632", "-0.2958"),
    (3, "0", "0.6752", "-0.0965"),
    (3, "1", "0.6604", "-0.2923"),
    (3, "2", "0.6348", "-0.4941"),
    (4, "0", "0.8673", "-0.0964"),
    (4, "1", "0.8857", "-0.2909"),
    (4, "2", "0.8345", "-0.4895"),
    (4, "3", "0.8064", "-0.6926"),
    (5, "0", "1.0585", "-0.096225"),
    (5, "1", "1.0585", "-0.28868"),
    (5, "2", "1.0585", "-0.48113"),
    (5, "3", "1.0585", "-0.67358"),
    (5, "4", "1.0585", "-0.86603"),
    (5, "5", "1.0585", "-1.0585"),
    (6, "0", "1.2509", "-0.096225"),
    (6, "1", "1.2509", "-0.28868"),
    (6, "2", "1.2509", "-0.48113"),
    (6, "3", "1.2509", "-0.67358"),
    (6, "4", "1.2509", "-0.86603"),
    (6, "5", "1.2509", "-1.0585"),
    (6, "6", "1.2509", "-1.2509"),
    (7, "0", "1.4434", "-0.096225"),
    (7, "1", "1.4434", "-0.28868"),
    (7, "2", "1.4434", "-0.48113"),
    (7, "3", "1.4434", "-0.67358"),
    (7, "4", "1.4434", "-0.86603"),
    (7, "5", "1.4434", "-1.0585"),
    (7, "6", "1.4434", "-1.2509"),
    (7, "7", "1.4434", "-1.4434"),
    (8, "0", "1.6358", "-0.096225"),
    (8, "1", "1.6358", "-0.28868"),
    (8, "2", "1.6358", "-0.48113"),
    (8, "3", "1.6358", "-0.67358"),
    (8, "4", "1.6358", "-0.86603"),
    (8, "5", "1.6358", "-1.0585"),
    (8, "6", "1.6358", "-1.2509"),
    (8, "7", "1.6358", "-1.4434"),
    (8, "8", "1.6358", "-1.6358"),
    (9, "0", "1.8283", "-0.096225"),
    (9, "1", "1.8283", "-0.28868"),
    (9, "2", "1.8283", "-0.48113"),
    (9, "3", "1.8283", "-0.67358"),
    (9, "4", "1.8283", "-0.86603"),
    (9, "5", "1.8283", "-1.0585"),
    (9, "6", "1.8283", "-1.2509"),
    (9, "7", "1.8283", "-1.4434"),
    (9, "8", "1.8283", "-1.6358"),
    (9, "9", "1.8283", "-1.8283"),
    (10, "0", "2.0207", "-0.096225"),
    (10, "1", "2.0207", "-0.28868"),
    (10, "2", "2.0207", "-0.48113"),
    (10, "3", "2.0207", "-0.67358"),
    (10, "4", "2.0207", "-0.86603"),
    (10, "5", "2.0207", "-1.0585"),
    (10, "6", "2.0207", "-1.2509"),
    (10, "7", "2.0207", "-1.4434"),
    (10, "8", "2.0207", "-1.6358"),
    (10, "9", "2.0207", "-1.8283"),
    (10, "10", "2.0207", "-2.0207"),
    (11, "0", "2.2132", "-0.096225"),
    (11, "1", "2.2132", "-0.28868"),
    (11, "2", "2.2132", "-0.48113"),
    (11, "3", "2.2132", "-0.67358"),
    (11, "4", "2.2132", "-0.86603"),
    (11, "5", "2.2132", "-1.0585"),
    (11, "6", "2.2132", "-1.2509"),
    (11, "7", "2.2132", "-1.4434"),
    (11, "8", "2.2132", "-1.6358"),
    (11, "9", "2.2132", "-1.8283"),
    (11, "10", "2.2132", "-2.0207"),
    (11, "11", "2.2132", "-2.2132"),
    (12, "0", "2.4056", "-0.096225"),
    (12, "1", "2.4056", "-0.28868"),
    (12, "2", "2.4056", "-0.48113"),
    (12, "3", "2.4056", "-0.67358"),
    (12, "4", "2.4056", "-0.86603"),
    (12, "5", "2.4056", "-1.0585"),
    (12, "6", "2.4056", "-1.2509"),
    (12, "7", "2.4056", "-1.4434"),
    (12, "8", "2.4056", "-1.6358"),
    (12, "9", "2.4056", "-1.8283"),
    (12, "10", "2.4056", "-2.0207"),
    (12, "11", "2.4056", "-2.2132"),
    (12, "12", "2.4056", "-2.4056"),
    (13, "0", "2.5981", "-0.096225"),
    (13, "1", "2.5981", "-0.28868"),
    (13, "2", "2.5981", "-0.48113"),
    (13, "3", "2.5981", "-0.67358"),
    (13, "4", "2.5981", "-0.86603"),
    (13, "5", "2.5981", "-1.0585"),
    (13, "6", "2.5981", "-1.2509"),
    (13, "7", "2.5981", "-1.4434"),
    (13, "8", "2.5981", "-1.6358"),
    (13, "9", "2.5981", "-1.8283"),
    (13, "10", "2.5981", "-2.0207"),
    (13, "11", "2.5981", "-2.2132"),
    (13, "12", "2.5981", "-2.4056"),
    (13, "13", "2.5981", "-2.5981"),
    (14, "0", "2.7905", "-0.096225"),
    (14, "1", "2.7905", "-0.28868"),
    (14, "2", "2.7905", "-0.48113"),
    (14, "3", "2.7905", "-0.67358"),
    (14, "4", "2.7905", "-0.86603"),
    (14, "5", "2.7905", "-1.0585"),
    (14, "6", "2.7905", "-1.2509"),
    (14, "7", "2.7905", "-1.4434"),
    (14, "8", "2.7905", "-1.6358"),
    (14, "9", "2.7905", "-1.8283"),
    (14, "10", "2.7905", "-2.0207"),
    (14, "11", "2.7905", "-2.2132"),
    (14, "12", "2.7905", "-2.4056"),
    (14, "13", "2.7905", "-2.5981"),
    (14, "14", "2.7905", "-2.7905"),
    (15, "0", "2.9830", "-0.096225"),
    (15, "1", "2.9830", "-0.28868"),
    (15, "2", "2.9830", "-0.48113"),
    (15, "3", "2.9830", "-0.67358"),
    (15, "4", "2.9830", "-0.86603"),
    (15, "5", "2.9830", "-1.0585"),
    (15, "6", "2.9830", "-1.2509"),
    (15, "7", "2.9830", "-1.4434"),
    (15, "8", "2.9830", "-1.6358"),
    (15, "9", "2.9830", "-1.8283"),
    (15, "10", "2.9830", "-2.0207"),
    (15, "11", "2.9830", "-2.2132"),
    (15, "12", "2.9830", "-2.4056"),
    (15, "13", "2.9830", "-2.5981"),
    (15, "14", "2.9830", "-2.7905"),
    (15, "15", "2.9830", "-2.9830"),
    (16, "0", "3.1754", "-0.096225"),
    (16, "1", "3.1754", "-0.28868"),
    (16, "2", "3.1754", "-0.48113"),
    (16, "3", "3.1754", "-0.67358"),
    (16, "4", "3.1754", "-0.86603"),
    (16, "5", "3.1754", "-1.0585"),
    (16, "6", "3.1754", "-1.2509"),
    (16, "7", "3.1754", "-1.4434"),
    (16, "8", "3.1754", "-1.6358"),
    (16, "9", "3.1754", "-1.8283"),
    (16, "10", "3.1754", "-2.0207"),
    (16, "11", "3.1754", "-2.2132"),
    (16, "12", "3.1754", "-2.4056"),
    (16, "13", "3.1754", "-2.5981"),
    (16, "14", "3.1754", "-2.7905"),
    (16, "15", "3.1754", "-2.9830"),
    (16, "16", "3.1754", "-3.1754"),
    (17, "0", "3.3679", "-0.096225"),
    (17, "1", "3.3679", "-0.28868"),
    (17, "2", "3.3679", "-0.48113"),
    (17, "3", "3.3679", "-0.67358"),
    (17, "4", "3.3679", "-0.86603"),
    (17, "5", "3.3679", "-1.0585"),
    (17, "6", "3.3679", "-1.2509"),
    (17, "7", "3.3679", "-1.4434"),
    (17, "8", "3.3679", "-1.6358"),
    (17, "9", "3.3679", "-1.8283"),
    (17, "10", "3.3679", "-2.0207"),
    (17, "11", "3.3679", "-2.2132"),
    (17, "12", "3.3679", "-2.4056"),
    (17, "13", "3.3679", "-2.5981"),
    (17, "14", "3.3679", "-2.7905"),
    (17, "15", "3.3679", "-2.9830"),
    (17, "16", "3.3679", "-3.1754"),
    (17, "17", "3.3679", "-3.3679"),
    (18, "0", "3.5603", "-0.096225"),
    (18, "1", "3.5603", "-0.28868"),
    (18, "2", "3.5603", "-0.48113"),
    (18, "3", "3.5603", "-0.67358"),
    (18, "4", "3.5603", "-0.86603"),
    (18, "5", "3.5603", "-1.0585"),
    (18, "6", "3.5603", "-1.2509"),
    (18, "7", "3.5603", "-1.4434"),
    (18, "8", "3.5603", "-1.6358"),
    (18, "9", "3.5603", "-1.8283"),
    (18, "10", "3.5603", "-2.0207"),
    (18, "11", "3.5603", "-2.2132"),
    (18, "12", "3.5603", "-2.4056"),
    (18, "13", "3.5603", "-2.5981"),
    (18, "14", "3.5603", "-2.7905"),
    (18, "15", "3.5603", "-2.9830"),
    (18, "16", "3.5603", "-3.1754"),
    (18, "17", "3.5603", "-3.3679"),
    (18, "18", "3.5603", "-3.5603"),
    (19, "0", "3.7528", "-0.096225"),
    (19, "1", "3.7528", "-0.28868"),
    (19, "2", "3.7528", "-0.48113"),
    (19, "3", "3.7528", "-0.67358"),
    (19, "4", "3.7528", "-0.86603"),
    (19, "5", "3.7528", "-1.0585"),
    (19, "6", "3.7528", "-1.2509"),
    (19, "7", "3.7528", "-1.4434"),
    (19, "8", "3.7528", "-1.6358"),
    (19, "9", "3.7528", "-1.8283"),
    (19, "10", "3.7528", "-2.0207"),
    (19, "11", "3.7528", "-2.2132"),
    (19, "12", "3.7528", "-2.4056"),
    (19, "13", "3.7528", "-2.5981"),
    (19, "14", "3.7528", "-2.7905"),
    (19, "15", "3.7528", "-2.9830"),
    (19, "16", "3.7528", "-3.1754"),
    (19, "17", "3.7528", "-3.3679"),
    (19, "18", "3.7528", "-3.5603"),
    (19, "19", "3.7528", "-3.7528"),
    (20, "0", "3.9452", "-0.096225"),
    (20, "1", "3.9452", "-0.28868"),
    (20, "2", "3.9452", "-0.48113"),
    (20, "3", "3.9452", "-0.67358"),
    (20, "4", "3.9452", "-0.86603"),
    (20, "5", "3.9452", "-1.0585"),
    (20, "6", "3.9452", "-1.2509"),
    (20, "7", "3.9452", "-1.4434"),
    (20, "8", "3.9452", "-1.6358"),
    (20, "9", "3.9452", "-1.8283"),
    (20, "10", "3.9452", "-2.0207"),
    (20, "11", "3.9452", "-2.2132"),
    (20, "12", "3.9452", "-2.4056"),
    (20, "13", "3.9452", "-2.5981"),
    (20, "14", "3.9452", "-2.7905"),
    (20, "15", "3.9452", "-2.9830"),
    (20, "16", "3.9452", "-3.1754"),
    (20, "17", "3.9452", "-3.3679"),
    (20, "18", "3.9452", "-3.5603"),
    (20, "19", "3.9452", "-3.7528"),
    (20, "20", "3.9452", "-3.9452"),
    (21, "0", "4.1377", "-0.096225"),
    (21, "1", "4.1377", "-0.28868"),
    (21, "2", "4.1377", "-0.48113"),
    (21, "3", "4.1377", "-0.67358"),
    (21, "4", "4.1377", "-0.86603"),
    (21, "5", "4.1377", "-1.0585"),
    (21, "6", "4.1377", "-1.2509"),
    (21, "7", "4.1377", "-1.4434"),
    (21, "8", "4.1377", "-1.6358"),
    (21, "9", "4.1377", "-1.8283"),
    (21, "10", "4.1377", "-2.0207"),
    (21, "11", "4.1377", "-2.2132"),
    (21, "12", "4.1377", "-2.4056"),
    (21, "13", "4.1377", "-2.5981"),
    (21, "14", "4.1377", "-2.7905"),
    (21, "15", "4.1377", "-2.9830"),
    (21, "16", "4.1377", "-3.1754"),
    (21, "17", "4.1377", "-3.3679"),
    (21, "18", "4.1377", "-3.5603"),
    (21, "19", "4.1377", "-3.7528"),
    (21, "20", "4.1377", "-3.9452"),
    (21, "21", "4.1377", "-4.1377"),
    (22, "0", "4.3301", "-0.096225"),
    (22, "1", "4.3301", "-0.28868"),
    (22, "2", "4.3301", "-0.48113"),
    (22, "3", "4.3301", "-0.67358"),
    (22, "4", "4.3301", "-0.86603"),
    (22, "5", "4.3301", "-1.0585"),
    (22, "6", "4.3301", "-1.2509"),
    (22, "7", "4.3301", "-1.4434"),
    (22, "8", "4.3301", "-1.6358"),
    (22, "9", "4.3301", "-1.8283"),
    (22, "10", "4.3301", "-2.0207"),
    (22, "11", "4.3301", "-2.2132"),
    (22, "12", "4.3301", "-2.4056"),
    (22, "13", "4.3301", "-2.5981"),
    (22, "14", "4.3301", "-2.7905"),
    (22, "15", "4.3301", "-2.9830"),
    (22, "16", "4.3301", "-3.1754"),
    (22, "17", "4.3301", "-3.3679"),
    (22, "18", "4.3301", "-3.5603"),
    (22, "19", "4.3301", "-3.7528"),
    (22, "20", "4.3301", "-3.9452"),
    (22, "21", "4.3301", "-4.1377"),
    (22, "22", "4.3301", "-4.3301"),
    (23, "0", "4.5226", "-0.096225"),
    (23, "1", "4.5226", "-0.28868"),
    (23, "2", "4.5226", "-0.48113"),
    (23, "3", "4.5226", "-0.67358"),
    (23, "4", "4.5226", "-0.86603"),
    (23, "5", "4.5226", "-1.0585"),
    (23, "6", "4.5226", "-1.2509"),
    (23, "7", "4.5226", "-1.4434"),
    (23, "8", "4.5226", "-1.6358"),
    (23, "9", "4.5226", "-1.8283"),
    (23, "10", "4.5226", "-2.0207"),
    (23, "11", "4.5226", "-2.2132"),
    (23, "12", "4.5226", "-2.4056"),
    (23, "13", "4.5226", "-2.5981"),
    (23, "14", "4.5226", "-2.7905"),
    (23, "15", "4.5226", "-2.9830"),
    (23, "16", "4.5226", "-3.1754"),
    (23, "17", "4.5226", "-3.3679"),
    (23, "18", "4.5226", "-3.5603"),
    (23, "19", "4.5226", "-3.7528"),
    (23, "20", "4.5226", "-3.9452"),
    (23, "21", "4.5226", "-4.1377"),
    (23, "22", "4.5226", "-4.3301"),
    (23, "23", "4.5226", "-4.5226"),
    (24, "0", "4.7150", "-0.096225"),
    (24, "1", "4.7150", "-0.28868"),
    (24, "2", "4.7150", "-0.48113"),
    (24, "3", "4.7150", "-0.67358"),
    (24, "4", "4.7150", "-0.86603"),
    (24, "5", "4.7150", "-1.0585"),
    (24, "6", "4.7150", "-1.2509"),
    (24, "7", "4.7150", "-1.4434"),
    (24, "8", "4.7150", "-1.6358"),
    (24, "9", "4.7150", "-1.8283"),
    (24, "10", "4.7150", "-2.0207"),
    (24, "11", "4.7150", "-2.2132"),
    (24, "12", "4.7150", "-2.4056"),
    (24, "13", "4.7150", "-2.5981"),
    (24, "14", "4.7150", "-2.7905"),
    (24, "15", "4.7150", "-2.9830"),
    (24, "16", "4.7150", "-3.1754"),
    (24, "17", "4.7150", "-3.3679"),
    (24, "18", "4.7150", "-3.5603"),
    (24, "19", "4.7150", "-3.7528"),
    (24, "20", "4.7150", "-3.9452"),
    (24, "21", "4.7150", "-4.1377"),
    (24, "22", "4.7150", "-4.3301"),
    (24, "23", "4.7150", "-4.5226"),
    (24, "24", "4.7150", "-4.7150"),
    (25, "0", "4.9075", "-0.096225"),
    (25, "1", "4.9075", "-0.28868"),
    (25, "2", "4.9075", "-0.48113"),
    (25, "3", "4.9075", "-0.67358"),
    (25, "4", "4.9075", "-0.86603"),
    (25, "5", "4.9075", "-1.0585"),
    (25, "6", "4.9075", "-1.2509"),
    (25, "7", "4.9075", "-1.4434"),
    (25, "8", "4.9075", "-1.6358"),
    (25, "9", "4.9075", "-1.8283"),
    (25, "10", "4.9075", "-2.0207"),
    (25, "11", "4.9075", "-2.2132"),
    (25, "12", "4.9075", "-2.4056"),
    (25, "13", "4.9075", "-2.5981"),
    (25, "14", "4.9075", "-2.7905"),
    (25, "15", "4.9075", "-2.9830"),
    (25, "16", "4.9075", "-3.1754"),
    (25, "17", "4.9075", "-3.3679"),
    (25, "18", "4.9075", "-3.5603"),
    (25, "19", "4.9075", "-3.7528"),
    (25, "20", "4.9075", "-3.9452"),
    (25, "21", "4.9075", "-4.1377"),
    (25, "22", "4.9075", "-4.3301"),
    (25, "23", "4.9075", "-4.5226"),
    (25, "24", "4.9075", "-4.7150"),
    (25, "25", "4.9075", "-4.9075"),
    (26, "0", "5.0999", "-0.096225"),
    (26, "1", "5.0999", "-0.28868"),
    (26, "2", "5.0999", "-0.48113"),
    (26, "3", "5.0999", "-0.67358"),
    (26, "4", "5.0999", "-0.86603"),
    (26, "5", "5.0999", "-1.0585"),
    (26, "6", "5.0999", "-1.2509"),
    (26, "7", "5.0999", "-1.4434"),
    (26, "8", "5.0999", "-1.6358"),
    (26, "9", "5.0999", "-1.8283"),
    (26, "10", "5.0999", "-2.0207"),
    (26, "11", "5.0999", "-2.2132"),
    (26, "12", "5.0999", "-2.4056"),
    (26, "13", "5.0999", "-2.5981"),
    (26, "14", "5.0999", "-2.7905"),
    (26, "15", "5.0999", "-2.9830"),
    (26, "16", "5.0999", "-3.1754"),
    (26, "17", "5.0999", "-3.3679"),
    (26, "18", "5.0999", "-3.5603"),
    (26, "19", "5.0999", "-3.7528"),
    (26, "20", "5.0999", "-3.9452"),
    (26, "21", "5.0999", "-4.1377"),
    (26, "22", "5.0999", "-4.3301"),
    (26, "23", "5.0999", "-4.5226"),
    (26, "24", "5.0999", "-4.7150"),
    (26, "25", "5.0999", "-4.9075"),
    (26, "26", "5.0999", "-5.0999"),
    (27, "0", "5.2924", "-0.096225"),
    (27, "1", "5.2924", "-0.28868"),
    (27, "2", "5.2924", "-0.48113"),
    (27, "3", "5.2924", "-0.67358"),
    (27, "4", "5.2924", "-0.86603"),
    (27, "5", "5.2924", "-1.0585"),
    (27, "6", "5.2924", "-1.2509"),
    (27, "7", "5.2924", "-1.4434"),
    (27, "8", "5.2924", "-1.6358"),
    (27, "9", "5.2924", "-1.8283"),
    (27, "10", "5.2924", "-2.0207"),
    (27, "11", "5.2924", "-2.2132"),
    (27, "12", "5.2924", "-2.4056"),
    (27, "13", "5.2924", "-2.5981"),
    (27, "14", "5.2924", "-2.7905"),
    (27, "15", "5.2924", "-2.9830"),
    (27, "16", "5.2924", "-3.1754"),
    (27, "17", "5.2924", "-3.3679"),
    (27, "18", "5.2924", "-3.5603"),
    (27, "19", "5.2924", "-3.7528"),
    (27, "20", "5.2924", "-3.9452"),
    (27, "21", "5.2924", "-4.1377"),
    (27, "22", "5.2924", "-4.3301"),
    (27, "23", "5.2924", "-4.5226"),
    (27, "24", "5.2924", "-4.7150"),
    (27, "25", "5.2924", "-4.9075"),
    (27, "26", "5.2924", "-5.0999"),
    (27, "27", "5.2924", "-5.2924"),
    (28, "0", "5.4848", "-0.096225"),
    (28, "1", "5.4848", "-0.28868"),
    (28, "2", "5.4848", "-0.48113"),
    (28, "3", "5.4848", "-0.67358"),
    (28, "4", "5.4848", "-0.86603"),
    (28, "5", "5.4848", "-1.0585"),
    (28, "6", "5.4848", "-1.2509"),
    (28, "7", "5.4848", "-1.4434"),
    (28, "8", "5.4848", "-1.6358"),
    (28, "9", "5.4848", "-1.8283"),
    (28, "10", "5.4848", "-2.0207"),
    (28, "11", "5.4848", "-2.2132"),
    (28, "12", "5.4848", "-2.4056"),
    (28, "13", "5.4848", "-2.5981"),
    (28, "14", "5.4848", "-2.7905"),
    (28, "15", "5.4848", "-2.9830"),
    (28, "16", "5.4848", "-3.1754"),
    (28, "17", "5.4848", "-3.3679"),
    (28, "18", "5.4848", "-3.5603"),
    (28, "19", "5.4848", "-3.7528"),
    (28, "20", "5.4848", "-3.9452"),
    (28, "21", "5.4848", "-4.1377"),
    (28, "22", "5.4848", "-4.3301"),
    (28, "23", "5.4848", "-4.5226"),
    (28, "24", "5.4848", "-4.7150"),
    (28, "25", "5.4848", "-4.9075"),
    (28, "26", "5.4848", "-5.0999"),
    (28, "27", "5.4848", "-5.2924"),
    (28, "28", "5.4848", "-5.4848"),
    (29, "0", "5.6773", "-0.096225"),
    (29, "1", "5.6773", "-0.28868"),
    (29, "2", "5.6773", "-0.48113"),
    (29, "3", "5.6773", "-0.67358"),
    (29, "4", "5.6773", "-0.86603"),
    (29, "5", "5.6773", "-1.0585"),
    (29, "6", "5.6773", "-1.2509"),
    (29, "7", "5.6773", "-1.4434"),
    (29, "8", "5.6773", "-1.6358"),
    (29, "9", "5.6773", "-1.8283"),
    (29, "10", "5.6773", "-2.0207"),
    (29, "11", "5.6773", "-2.2132"),
    (29, "12", "5.6773", "-2.4056"),
    (29, "13", "5.6773", "-2.5981"),
    (29, "14", "5.6773", "-2.7905"),
    (29, "15", "5.6773", "-2.9830"),
    (29, "16", "5.6773", "-3.1754"),
    (29, "17", "5.6773", "-3.3679"),
    (29, "18", "5.6773", "-3.5603"),
    (29, "19", "5.6773", "-3.7528"),
    (29, "20", "5.6773", "-3.9452"),
    (29, "21", "5.6773", "-4.1377"),
    (29, "22", "5.6773", "-4.3301"),
    (29, "23", "5.6773", "-4.5226"),
    (29, "24", "5.6773", "-4.7150"),
    (29, "25", "5.6773", "-4.9075"),
    (29, "26", "5.6773", "-5.0999"),
    (29, "27", "5.6773", "-5.2924"),
    (29, "28", "5.6773", "-5.4848"),
    (29, "29", "5.6773", "-5.6773"),
    (30, "0", "5.8697", "-0.096225"),
    (30, "1", "5.8697", "-0.28868"),
    (30, "2", "5.8697", "-0.48113"),
    (30, "3", "5.8697", "-0.67358"),
    (30, "4", "5.8697", "-0.86603"),
    (30, "5", "5.8697", "-1.0585"),
    (30, "6", "5.8697", "-1.2509"),
    (30, "7", "5.8697", "-1.4434"),
    (30, "8", "5.8697", "-1.6358"),
    (30, "9", "5.8697", "-1.8283"),
    (30, "10", "5.8697", "-2.0207"),
    (30, "11", "5.8697", "-2.2132"),
    (30, "12", "5.8697", "-2.4056"),
    (30, "13", "5.8697", "-2.5981"),
    (30, "14", "5.8697", "-2.7905"),
    (30, "15", "5.8697", "-2.9830"),
    (30, "16", "5.8697", "-3.1754"),
    (30, "17", "5.8697", "-3.3679"),
    (30, "18", "5.8697", "-3.5603"),
    (30, "19", "5.8697", "-3.7528"),
    (30, "20", "5.8697", "-3.9452"),
    (30, "21", "5.8697", "-4.1377"),
    (30, "22", "5.8697", "-4.3301"),
    (30, "23", "5.8697", "-4.5226"),
    (30, "24", "5.8697", "-4.7150"),
    (30, "25", "5.8697", "-4.9075"),
    (30, "26", "5.8697", "-5.0999"),
    (30, "27", "5.8697", "-5.2924"),
    (30, "28", "5.8697", "-5.4848"),
    (30, "29", "5.8697", "-5.6773"),
    (30, "30", "5.8697", "-5.8697"),
    (31,"0","6.0622","-0.096225"),
    (31,"1","6.0622","-0.28868"),
    (31,"2","6.0622","-0.48113"),
    (31,"3","6.0622","-0.67358"),
    (31,"4","6.0622","-0.86603"),
    (31,"5","6.0622","-1.0585"),
    (31,"6","6.0622","-1.2509"),
    (31,"7","6.0622","-1.4434"),
    (31,"8","6.0622","-1.6358"),
    (31,"9","6.0622","-1.8283"),
    (31,"10","6.0622","-2.0207"),
    (31,"11","6.0622","-2.2132"),
    (31,"12","6.0622","-2.4056"),
    (31,"13","6.0622","-2.5981"),
    (31,"14","6.0622","-2.7905"),
    (31,"15","6.0622","-2.9830"),
    (31,"16","6.0622","-3.1754"),
    (31,"17","6.0622","-3.3679"),
    (31,"18","6.0622","-3.5603"),
    (31,"19","6.0622","-3.7528"),
    (31,"20","6.0622","-3.9452"),
    (31,"21","6.0622","-4.1377"),
    (31,"22","6.0622","-4.3301"),
    (31,"23","6.0622","-4.5226"),
    (31,"24","6.0622","-4.7150"),
    (31,"25","6.0622","-4.9075"),
    (31,"26","6.0622","-5.0999"),
    (31,"27","6.0622","-5.2924"),
    (31,"28","6.0622","-5.4848"),
    (31,"29","6.0622","-5.6773"),
    (31,"30","6.0622","-5.8697"),
    (31,"31","6.0622","-6.0622"),
    (31,"32","6.0622","-6.2546"),
    (31,"33","6.0622","-6.4471"),
    (31,"34","6.0622","-6.6395"),
    (31,"35","6.0622","-6.8320"),
    (31,"36","6.0622","-7.0244"),
    (31,"37","6.0622","-7.2169"),
    (31,"38","6.0622","-7.4093"),
    (31,"39","6.0622","-7.6018"),
    (32,"0","6.2546","-0.096225"),
    (32,"1","6.2546","-0.28868"),
    (32,"2","6.2546","-0.48113"),
    (32,"3","6.2546","-0.67358"),
    (32,"4","6.2546","-0.86603"),
    (32,"5","6.2546","-1.0585"),
    (32,"6","6.2546","-1.2509"),
    (32,"7","6.2546","-1.4434"),
    (32,"8","6.2546","-1.6358"),
    (32,"9","6.2546","-1.8283"),
    (32,"10","6.2546","-2.0207"),
    (32,"11","6.2546","-2.2132"),
    (32,"12","6.2546","-2.4056"),
    (32,"13","6.2546","-2.5981"),
    (32,"14","6.2546","-2.7905"),
    (32,"15","6.2546","-2.9830"),
    (32,"16","6.2546","-3.1754"),
    (32,"17","6.2546","-3.3679"),
    (32,"18","6.2546","-3.5603"),
    (32,"19","6.2546","-3.7528"),
    (32,"20","6.2546","-3.9452"),
    (32,"21","6.2546","-4.1377"),
    (32,"22","6.2546","-4.3301"),
    (32,"23","6.2546","-4.5226"),
    (32,"24","6.2546","-4.7150"),
    (32,"25","6.2546","-4.9075"),
    (32,"26","6.2546","-5.0999"),
    (32,"27","6.2546","-5.2924"),
    (32,"28","6.2546","-5.4848"),
    (32,"29","6.2546","-5.6773"),
    (32,"30","6.2546","-5.8697"),
    (32,"31","6.2546","-6.0622"),
    (32,"32","6.2546","-6.2546"),
    (32,"33","6.2546","-6.4471"),
    (32,"34","6.2546","-6.6395"),
    (32,"35","6.2546","-6.8320"),
    (32,"36","6.2546","-7.0244"),
    (32,"37","6.2546","-7.2169"),
    (32,"38","6.2546","-7.4093"),
    (32,"39","6.2546","-7.6018"),
    (33,"0","6.4471","-0.096225"),
    (33,"1","6.4471","-0.28868"),
    (33,"2","6.4471","-0.48113"),
    (33,"3","6.4471","-0.67358"),
    (33,"4","6.4471","-0.86603"),
    (33,"5","6.4471","-1.0585"),
    (33,"6","6.4471","-1.2509"),
    (33,"7","6.4471","-1.4434"),
    (33,"8","6.4471","-1.6358"),
    (33,"9","6.4471","-1.8283"),
    (33,"10","6.4471","-2.0207"),
    (33,"11","6.4471","-2.2132"),
    (33,"12","6.4471","-2.4056"),
    (33,"13","6.4471","-2.5981"),
    (33,"14","6.4471","-2.7905"),
    (33,"15","6.4471","-2.9830"),
    (33,"16","6.4471","-3.1754"),
    (33,"17","6.4471","-3.3679"),
    (33,"18","6.4471","-3.5603"),
    (33,"19","6.4471","-3.7528"),
    (33,"20","6.4471","-3.9452"),
    (33,"21","6.4471","-4.1377"),
    (33,"22","6.4471","-4.3301"),
    (33,"23","6.4471","-4.5226"),
    (33,"24","6.4471","-4.7150"),
    (33,"25","6.4471","-4.9075"),
    (33,"26","6.4471","-5.0999"),
    (33,"27","6.4471","-5.2924"),
    (33,"28","6.4471","-5.4848"),
    (33,"29","6.4471","-5.6773"),
    (33,"30","6.4471","-5.8697"),
    (33,"31","6.4471","-6.0622"),
    (33,"32","6.4471","-6.2546"),
    (33,"33","6.4471","-6.4471"),
    (33,"34","6.4471","-6.6395"),
    (33,"35","6.4471","-6.8320"),
    (33,"36","6.4471","-7.0244"),
    (33,"37","6.4471","-7.2169"),
    (33,"38","6.4471","-7.4093"),
    (33,"39","6.4471","-7.6018"),
    (34,"0","6.6395","-0.096225"),
    (34,"1","6.6395","-0.28868"),
    (34,"2","6.6395","-0.48113"),
    (34,"3","6.6395","-0.67358"),
    (34,"4","6.6395","-0.86603"),
    (34,"5","6.6395","-1.0585"),
    (34,"6","6.6395","-1.2509"),
    (34,"7","6.6395","-1.4434"),
    (34,"8","6.6395","-1.6358"),
    (34,"9","6.6395","-1.8283"),
    (34,"10","6.6395","-2.0207"),
    (34,"11","6.6395","-2.2132"),
    (34,"12","6.6395","-2.4056"),
    (34,"13","6.6395","-2.5981"),
    (34,"14","6.6395","-2.7905"),
    (34,"15","6.6395","-2.9830"),
    (34,"16","6.6395","-3.1754"),
    (34,"17","6.6395","-3.3679"),
    (34,"18","6.6395","-3.5603"),
    (34,"19","6.6395","-3.7528"),
    (34,"20","6.6395","-3.9452"),
    (34,"21","6.6395","-4.1377"),
    (34,"22","6.6395","-4.3301"),
    (34,"23","6.6395","-4.5226"),
    (34,"24","6.6395","-4.7150"),
    (34,"25","6.6395","-4.9075"),
    (34,"26","6.6395","-5.0999"),
    (34,"27","6.6395","-5.2924"),
    (34,"28","6.6395","-5.4848"),
    (34,"29","6.6395","-5.6773"),
    (34,"30","6.6395","-5.8697"),
    (34,"31","6.6395","-6.0622"),
    (34,"32","6.6395","-6.2546"),
    (34,"33","6.6395","-6.4471"),
    (34,"34","6.6395","-6.6395"),
    (34,"35","6.6395","-6.8320"),
    (34,"36","6.6395","-7.0244"),
    (34,"37","6.6395","-7.2169"),
    (34,"38","6.6395","-7.4093"),
    (34,"39","6.6395","-7.6018"),
    (35,"0","6.8320","-0.096225"),
    (35,"1","6.8320","-0.28868"),
    (35,"2","6.8320","-0.48113"),
    (35,"3","6.8320","-0.67358"),
    (35,"4","6.8320","-0.86603"),
    (35,"5","6.8320","-1.0585"),
    (35,"6","6.8320","-1.2509"),
    (35,"7","6.8320","-1.4434"),
    (35,"8","6.8320","-1.6358"),
    (35,"9","6.8320","-1.8283"),
    (35,"10","6.8320","-2.0207"),
    (35,"11","6.8320","-2.2132"),
    (35,"12","6.8320","-2.4056"),
    (35,"13","6.8320","-2.5981"),
    (35,"14","6.8320","-2.7905"),
    (35,"15","6.8320","-2.9830"),
    (35,"16","6.8320","-3.1754"),
    (35,"17","6.8320","-3.3679"),
    (35,"18","6.8320","-3.5603"),
    (35,"19","6.8320","-3.7528"),
    (35,"20","6.8320","-3.9452"),
    (35,"21","6.8320","-4.1377"),
    (35,"22","6.8320","-4.3301"),
    (35,"23","6.8320","-4.5226"),
    (35,"24","6.8320","-4.7150"),
    (35,"25","6.8320","-4.9075"),
    (35,"26","6.8320","-5.0999"),
    (35,"27","6.8320","-5.2924"),
    (35,"28","6.8320","-5.4848"),
    (35,"29","6.8320","-5.6773"),
    (35,"30","6.8320","-5.8697"),
    (35,"31","6.8320","-6.0622"),
    (35,"32","6.8320","-6.2546"),
    (35,"33","6.8320","-6.4471"),
    (35,"34","6.8320","-6.6395"),
    (35,"35","6.8320","-6.8320"),
    (35,"36","6.8320","-7.0244"),
    (35,"37","6.8320","-7.2169"),
    (35,"38","6.8320","-7.4093"),
    (35,"39","6.8320","-7.6018"),
    (36,"0","7.0244","-0.096225"),
    (36,"1","7.0244","-0.28868"),
    (36,"2","7.0244","-0.48113"),
    (36,"3","7.0244","-0.67358"),
    (36,"4","7.0244","-0.86603"),
    (36,"5","7.0244","-1.0585"),
    (36,"6","7.0244","-1.2509"),
    (36,"7","7.0244","-1.4434"),
    (36,"8","7.0244","-1.6358"),
    (36,"9","7.0244","-1.8283"),
    (36,"10","7.0244","-2.0207"),
    (36,"11","7.0244","-2.2132"),
    (36,"12","7.0244","-2.4056"),
    (36,"13","7.0244","-2.5981"),
    (36,"14","7.0244","-2.7905"),
    (36,"15","7.0244","-2.9830"),
    (36,"16","7.0244","-3.1754"),
    (36,"17","7.0244","-3.3679"),
    (36,"18","7.0244","-3.5603"),
    (36,"19","7.0244","-3.7528"),
    (36,"20","7.0244","-3.9452"),
    (36,"21","7.0244","-4.1377"),
    (36,"22","7.0244","-4.3301"),
    (36,"23","7.0244","-4.5226"),
    (36,"24","7.0244","-4.7150"),
    (36,"25","7.0244","-4.9075"),
    (36,"26","7.0244","-5.0999"),
    (36,"27","7.0244","-5.2924"),
    (36,"28","7.0244","-5.4848"),
    (36,"29","7.0244","-5.6773"),
    (36,"30","7.0244","-5.8697"),
    (36,"31","7.0244","-6.0622"),
    (36,"32","7.0244","-6.2546"),
    (36,"33","7.0244","-6.4471"),
    (36,"34","7.0244","-6.6395"),
    (36,"35","7.0244","-6.8320"),
    (36,"36","7.0244","-7.0244"),
    (36,"37","7.0244","-7.2169"),
    (36,"38","7.0244","-7.4093"),
    (36,"39","7.0244","-7.6018"),
    (37,"0","7.2169","-0.096225"),
    (37,"1","7.2169","-0.28868"),
    (37,"2","7.2169","-0.48113"),
    (37,"3","7.2169","-0.67358"),
    (37,"4","7.2169","-0.86603"),
    (37,"5","7.2169","-1.0585"),
    (37,"6","7.2169","-1.2509"),
    (37,"7","7.2169","-1.4434"),
    (37,"8","7.2169","-1.6358"),
    (37,"9","7.2169","-1.8283"),
    (37,"10","7.2169","-2.0207"),
    (37,"11","7.2169","-2.2132"),
    (37,"12","7.2169","-2.4056"),
    (37,"13","7.2169","-2.5981"),
    (37,"14","7.2169","-2.7905"),
    (37,"15","7.2169","-2.9830"),
    (37,"16","7.2169","-3.1754"),
    (37,"17","7.2169","-3.3679"),
    (37,"18","7.2169","-3.5603"),
    (37,"19","7.2169","-3.7528"),
    (37,"20","7.2169","-3.9452"),
    (37,"21","7.2169","-4.1377"),
    (37,"22","7.2169","-4.3301"),
    (37,"23","7.2169","-4.5226"),
    (37,"24","7.2169","-4.7150"),
    (37,"25","7.2169","-4.9075"),
    (37,"26","7.2169","-5.0999"),
    (37,"27","7.2169","-5.2924"),
    (37,"28","7.2169","-5.4848"),
    (37,"29","7.2169","-5.6773"),
    (37,"30","7.2169","-5.8697"),
    (37,"31","7.2169","-6.0622"),
    (37,"32","7.2169","-6.2546"),
    (37,"33","7.2169","-6.4471"),
    (37,"34","7.2169","-6.6395"),
    (37,"35","7.2169","-6.8320"),
    (37,"36","7.2169","-7.0244"),
    (37,"37","7.2169","-7.2169"),
    (37,"38","7.2169","-7.4093"),
    (37,"39","7.2169","-7.6018"),
    (38,"0","7.4093","-0.096225"),
    (38,"1","7.4093","-0.28868"),
    (38,"2","7.4093","-0.48113"),
    (38,"3","7.4093","-0.67358"),
    (38,"4","7.4093","-0.86603"),
    (38,"5","7.4093","-1.0585"),
    (38,"6","7.4093","-1.2509"),
    (38,"7","7.4093","-1.4434"),
    (38,"8","7.4093","-1.6358"),
    (38,"9","7.4093","-1.8283"),
    (38,"10","7.4093","-2.0207"),
    (38,"11","7.4093","-2.2132"),
    (38,"12","7.4093","-2.4056"),
    (38,"13","7.4093","-2.5981"),
    (38,"14","7.4093","-2.7905"),
    (38,"15","7.4093","-2.9830"),
    (38,"16","7.4093","-3.1754"),
    (38,"17","7.4093","-3.3679"),
    (38,"18","7.4093","-3.5603"),
    (38,"19","7.4093","-3.7528"),
    (38,"20","7.4093","-3.9452"),
    (38,"21","7.4093","-4.1377"),
    (38,"22","7.4093","-4.3301"),
    (38,"23","7.4093","-4.5226"),
    (38,"24","7.4093","-4.7150"),
    (38,"25","7.4093","-4.9075"),
    (38,"26","7.4093","-5.0999"),
    (38,"27","7.4093","-5.2924"),
    (38,"28","7.4093","-5.4848"),
    (38,"29","7.4093","-5.6773"),
    (38,"30","7.4093","-5.8697"),
    (38,"31","7.4093","-6.0622"),
    (38,"32","7.4093","-6.2546"),
    (38,"33","7.4093","-6.4471"),
    (38,"34","7.4093","-6.6395"),
    (38,"35","7.4093","-6.8320"),
    (38,"36","7.4093","-7.0244"),
    (38,"37","7.4093","-7.2169"),
    (38,"38","7.4093","-7.4093"),
    (38,"39","7.4093","-7.6018"),
    (39,"0","7.6018","-0.096225"),
    (39,"1","7.6018","-0.28868"),
    (39,"2","7.6018","-0.48113"),
    (39,"3","7.6018","-0.67358"),
    (39,"4","7.6018","-0.86603"),
    (39,"5","7.6018","-1.0585"),
    (39,"6","7.6018","-1.2509"),
    (39,"7","7.6018","-1.4434"),
    (39,"8","7.6018","-1.6358"),
    (39,"9","7.6018","-1.8283"),
    (39,"10","7.6018","-2.0207"),
    (39,"11","7.6018","-2.2132"),
    (39,"12","7.6018","-2.4056"),
    (39,"13","7.6018","-2.5981"),
    (39,"14","7.6018","-2.7905"),
    (39,"15","7.6018","-2.9830"),
    (39,"16","7.6018","-3.1754"),
    (39,"17","7.6018","-3.3679"),
    (39,"18","7.6018","-3.5603"),
    (39,"19","7.6018","-3.7528"),
    (39,"20","7.6018","-3.9452"),
    (39,"21","7.6018","-4.1377"),
    (39,"22","7.6018","-4.3301"),
    (39,"23","7.6018","-4.5226"),
    (39,"24","7.6018","-4.7150"),
    (39,"25","7.6018","-4.9075"),
    (39,"26","7.6018","-5.0999"),
    (39,"27","7.6018","-5.2924"),
    (39,"28","7.6018","-5.4848"),
    (39,"29","7.6018","-5.6773"),
    (39,"30","7.6018","-5.8697"),
    (39,"31","7.6018","-6.0622"),
    (39,"32","7.6018","-6.2546"),
    (39,"33","7.6018","-6.4471"),
    (39,"34","7.6018","-6.6395"),
    (39,"35","7.6018","-6.8320"),
    (39,"36","7.6018","-7.0244"),
    (39,"37","7.6018","-7.2169"),
    (39,"38","7.6018","-7.4093"),
    (39,"39","7.6018","-7.6018"),
    (40,"0","7.7942","-0.096225"),
    (40,"1","7.7942","-0.28868"),
    (40,"2","7.7942","-0.48113"),
    (40,"3","7.7942","-0.67358"),
    (40,"4","7.7942","-0.86603"),
    (40,"5","7.7942","-1.0585"),
    (40,"6","7.7942","-1.2509"),
    (40,"7","7.7942","-1.4434"),
    (40,"8","7.7942","-1.6358"),
    (40,"9","7.7942","-1.8283"),
    (40,"10","7.7942","-2.0207"),
    (40,"11","7.7942","-2.2132"),
    (40,"12","7.7942","-2.4056"),
    (40,"13","7.7942","-2.5981"),
    (40,"14","7.7942","-2.7905"),
    (40,"15","7.7942","-2.9830"),
    (40,"16","7.7942","-3.1754"),
    (40,"17","7.7942","-3.3679"),
    (40,"18","7.7942","-3.5603"),
    (40,"19","7.7942","-3.7528"),
    (40,"20","7.7942","-3.9452"),
    (40,"21","7.7942","-4.1377"),
    (40,"22","7.7942","-4.3301"),
    (40,"23","7.7942","-4.5226"),
    (40,"24","7.7942","-4.7150"),
    (40,"25","7.7942","-4.9075"),
    (40,"26","7.7942","-5.0999"),
    (40,"27","7.7942","-5.2924"),
    (40,"28","7.7942","-5.4848"),
    (40,"29","7.7942","-5.6773"),
    (40,"30","7.7942","-5.8697"),
    (40,"31","7.7942","-6.0622"),
    (40,"32","7.7942","-6.2546"),
    (40,"33","7.7942","-6.4471"),
    (40,"34","7.7942","-6.6395"),
    (40,"35","7.7942","-6.8320"),
    (40,"36","7.7942","-7.0244"),
    (40,"37","7.7942","-7.2169"),
    (40,"38","7.7942","-7.4093"),
    (40,"39","7.7942","-7.6018")
]