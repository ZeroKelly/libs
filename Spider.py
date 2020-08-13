def compose_header(head):
    DIC = dict()
    for i in head.splitlines():
        tmp = i.split(': ')
        DIC[tmp[0]] = tmp[1]
    return DIC
