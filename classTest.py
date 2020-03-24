class BaseClass(object):
    
    def __init__(self,a,b):
        self.attribute1 = a
        self.attribute2 = b

    def showBaseAttributes(self):
        print("Attribute1 = {}, Attribute2 = {}".format(self.attribute1, 
                                                        self.attribute2))

class SubClass(BaseClass):

    def __init__(self,a=1,b=10,c=100):
        self.attribute3 = c
        super().__init__(a,b)

    def showSubAttributes(self):
        print("Attribute3 = {}".format(self.attribute3))


def main():
    a,b,c = 0.0,1.0,2.0
    baseClass = BaseClass(a,b)
    baseClass.showBaseAttributes()

    subClass = SubClass()
    subClass.showSubAttributes()

if __name__ == '__main__':
    main()