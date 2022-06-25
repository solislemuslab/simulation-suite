#ifndef MS_Events_hpp
#define MS_Events_hpp

#include <iostream>
#include <string>

enum MSEventType { split, join };

class MSEvent {
    public:
        MSEventType getEventType(void) { return eventType; }
        
        double getTime(void) { return time; }
        void setTime(double t) { time = t; }

        virtual std::string toString(void) { }
        void print(void) { std::cout << toString() << std::endl; }
        
    protected:
        MSEventType eventType;
        double time;
};

class MSJoinEvent:public MSEvent {
    public:
        MSJoinEvent(double t, int min, int maj) {
            time = t;
            majorTaxa = maj;
            minorTaxa = min;
        }
        MSJoinEvent(double t, std::string min, std::string maj) {
            time = t;
            majorTaxa = stoi(maj);
            minorTaxa = stoi(min);
        }
        MSJoinEvent(double t, int min, std::string maj) {
            time = t;
            majorTaxa = stoi(maj);
            minorTaxa = min;
        }

        int getMajorTaxa(void) { return majorTaxa; }
        void setMajorTaxa(int m) { majorTaxa = m; }

        int getMinorTaxa(void) { return minorTaxa; }
        void setMinorTaxa(int m) { minorTaxa = m; }

        std::string toString(void) { return "-ej " + std::to_string(time) + " " + std::to_string(minorTaxa) + " " + std::to_string(majorTaxa); }

    protected:
        int majorTaxa;
        int minorTaxa;
        MSEventType eventType = join;
};

class MSSplitEvent:public MSEvent {
    public:
        MSSplitEvent(double t, int tax, double gam) {
            time = t;
            taxa = tax;
            gamma = gam;
        }
        MSSplitEvent(double t, std::string tax, double gam) {
            time = t;
            taxa = stoi(tax);
            gamma = gam;
        }

        int getTaxa(void) { return taxa; }
        void setTaxa(int t) { taxa = t; }

        double getGamma(void) { return gamma; }
        void setGamma(double g) { gamma = g; }

        std::string toString(void) { return "-es " + std::to_string(time) + " " + std::to_string(taxa) + " " + std::to_string(gamma); }

    protected:
        int taxa;
        double gamma;
        MSEventType eventType = split;
};

#endif