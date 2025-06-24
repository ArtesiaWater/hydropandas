import requests


def test_gmn():
    
    bro_id = "GMN000000000001"
    url = f"https://publiek.broservices.nl/gm/gmn/v1/objects/{bro_id}"
    r = requests.get(url)

    r.raise_for_status()


# def test_gmw_characteristics()


#     url = "https://publiek.broservices.nl/gm/gmw/v1/characteristics/searches?"
#     r = requests.post(url, json=data)

#     r.raise_for_status()


def test_gmw():
    bro_id = "GMW000000036287"
    url = f"https://publiek.broservices.nl/gm/gmw/v1/objects/{bro_id}"
    r = requests.get(url)

    r.raise_for_status()


def test_gmw_relations():

    bro_id = "GMW000000036287"
    url = f"https://publiek.broservices.nl/gm/v1/gmw-relations/{bro_id}"
    r = requests.get(url)

    r.raise_for_status()

def test_gld():

    bro_id = "GLD000000008061"
    url = f"https://publiek.broservices.nl/gm/gld/v1/objects/{bro_id}"
    r = requests.get(url)

    r.raise_for_status()





