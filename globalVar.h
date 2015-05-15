#ifndef GLOBALVAR_H
#define GLOBALVAR_H

#include <string>
#include <iostream>

//Base class
class IGlobalVar
{
public:
	virtual std::string StringID() = 0; 
	virtual bool IsType(std::string) = 0;

	void Added(bool boolValue) {
		added = boolValue;
	}

	bool IsAdded() {
		return added;
	}

	void SetStringID(std::string n) {
		stringID = n;
	};



protected:
	std::string stringID; //name of variable
	std::string type; //variable type
	bool added = false;
};

template <class T>
class GlobalVar : public IGlobalVar
{
private:
	void SetType() {
		if (std::is_same<T, int>::value) type = "int";
		else if (std::is_same<T, double>::value) type = "double";
		else if (std::is_same<T, bool>::value) type = "bool";
	}
public:
	GlobalVar() { //Constructor calls function to define class type;
		SetType();
	}

	virtual std::string StringID() {
		return stringID;
	}
	virtual bool IsType(std::string typeCheck) {
		if (type == typeCheck) return true;
		else return false;
	};

	void SetValue(T n) {
		val = n;
	};

	T Value(void) {
		return val;
	};

protected:
	T val;
};

#endif