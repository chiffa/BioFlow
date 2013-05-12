'''
Created on Mar 25, 2013

@author: Andrei Kucharavy, achiffa @ github

This small example demonstrates why Python is NOT suited for Object Oriented Programming
due to it's random switching between referencing objects and new object creation
'''


#We consider names are UIDs
Parents2Children={'George':['Chloe','Baptiste'],'Alexandre':['Joseph','Melissa','Jonny'],'Serge':['Josphine']}


class Person:
    name='null'
    children=[]
      
    def __init__(self,name,children=[]):
        self.name=name
        self.children=children
    
    def addChild(self,childObject):
        self.children.append(childObject)
    
    def __repr__(self):
        return "<Person('%s','%s')>" % (self.name,
                                        len(self.children))

class Society:
    Name2Person={}
    
    def __init__(self,P2Ch):
        for Parent in Parents2Children.keys():
            self.Name2Person[Parent]=Person(Parent)
        for Children in Parents2Children.values():
            for Child in Children:
                self.Name2Person[Child]=Person(Child)


#<====================>
Soc=Society(Parents2Children)

print Soc.Name2Person

for ParentName in Parents2Children.keys():
    for ChildName in Parents2Children[ParentName]:
        Soc.Name2Person[ParentName].addChild(Soc.Name2Person[ChildName])

print Soc.Name2Person


'''
Expected Output:
{'Baptiste': <Person('Baptiste','0')>, 'Chloe': <Person('Chloe','0')>, 'Josphine': <Person('Josphine','0')>, 'Alexandre': <Person('Alexandre','0')>, 'Joseph': <Person('Joseph','0')>, 'Melissa': <Person('Melissa','0')>, 'Jonny': <Person('Jonny','0')>, 'George': <Person('George','0')>, 'Serge': <Person('Serge','0')>}
{'Baptiste': <Person('Baptiste','0')>, 'Chloe': <Person('Chloe','0')>, 'Josphine': <Person('Josphine','0')>, 'Alexandre': <Person('Alexandre','3')>, 'Joseph': <Person('Joseph','0')>, 'Melissa': <Person('Melissa','0')>, 'Jonny': <Person('Jonny','0')>, 'George': <Person('George','2')>, 'Serge': <Person('Serge','1')>}

Real Output:
{'Baptiste': <Person('Baptiste','0')>, 'Chloe': <Person('Chloe','0')>, 'Josphine': <Person('Josphine','0')>, 'Alexandre': <Person('Alexandre','0')>, 'Joseph': <Person('Joseph','0')>, 'Melissa': <Person('Melissa','0')>, 'Jonny': <Person('Jonny','0')>, 'George': <Person('George','0')>, 'Serge': <Person('Serge','0')>}
{'Baptiste': <Person('Baptiste','6')>, 'Chloe': <Person('Chloe','6')>, 'Josphine': <Person('Josphine','6')>, 'Alexandre': <Person('Alexandre','6')>, 'Joseph': <Person('Joseph','6')>, 'Melissa': <Person('Melissa','6')>, 'Jonny': <Person('Jonny','6')>, 'George': <Person('George','6')>, 'Serge': <Person('Serge','6')>}


apparently, for Python2.7 all the [] in __init__ sequence of a class 
are deemed as referring to the same object.

Mind = 'Blown'

Apparently, __init__ instance's parameters are created once and for all
and are not called again each time this function is called, as simple logic would suggest it.

In order to make it work, you have to shallow copy the void list each time it you call the
'default' values in __init__ method. for instance here, line 21 should be changed as following:

self.children=children => self.children=children[:]
''' 