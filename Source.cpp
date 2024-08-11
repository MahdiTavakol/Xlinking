/*
* A complicated program to add AGE into the Coarse grained collagen
* CG files.
* I do not guarantee that it works for your purpose!!!
* mahditavakol90@gmail.com
*/


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <random>
using namespace std;

enum Xlink {AGE,ENZ,BTH};
struct Coor
{
	float x;
	float y;
};
struct Bond
{
	int first;
	int second;
	int type; // 2 for AGE and 3 for ENZ
};
struct Angle
{
	int first;
	int second;
	int third;
	int type;
};
struct Atom
{
	int id;
	float x;
	float y;
	float z;
};
vector<Atom> atoms;

void writeAGETCLContact(const string& tclFileName, const string& trajectoryFileName,
	const string& contactFileName, const float& cutoff, const float& cutoffLo);
void writeENZTCLContact(const string& tclFileName, const string& dataFileName,
	const string& contactFileName, const float& cutoff, const float& cutoffLo);
void XLinker(const float& health, const Xlink& xlink,
	const string& dataFileName, const string& contactFileName,
	const string& contactFile2Name, const int&  initialBonds,
	int& numberOfMolecules, const int& molsInRadius,
	const int& beadsPerMolecule, vector<Bond>& AGEbonds,
	vector<Bond>& ENZbonds, vector<Angle>& AGEangles,
	vector<Angle>& ENZangles, const Coor& centralAxis,
	const float& radius, const float& AGEradius,
	ofstream& log, const bool& mode1 ,const bool& AGEInternal);
void writeTeminus(const string& dataFileName, int& numberOfMolecules);
void writeTeminus(const string& trajectoryFileName,
	const int& numberOfCollagenAtoms, const int& beadsPerMolecule);
void XYZWriter(const string& xyzFileName, vector<float> coors, const int& atomType);
int addAGE(const string& dataFileName, const string& contactFileName,
	const int& XTarget, int& numberOfMolecules, const int&  beadsPerMolecule,
	vector<Bond>& bonds, vector<Angle>& angles, const Coor& centralAxis,
	const float& radius, const float& AGEradius, const bool& AGEInternal);
int addENZ(const string& dataFileName, const string& contactFileName,
	const int& XTarget, int& numberOfMolecules, const int&  beadsPerMolecule,
	vector<Bond>& bonds, vector<Angle>& angles,
	const Xlink& xlink);
int measureBeadsPerMolecule(const string& dataFileName, int& numberOfCollagenAtoms, int& numberOfMolecules);
int measureMolPerRadius(const string& dataFileName, const float& AGEradius, const Coor& centralAxis);
void mergeBonds(vector<Bond>& bonds, const vector<Bond>& AGEbonds, const vector<Bond>& ENZbonds, const bool debug);
void mergeAngles(vector<Angle>& angles, const vector<Angle>& AGEangles, const vector<Angle>& ENZangles);
void writeDataFile(const string& dataFileName, const string& outputDataFileName,
	const vector<Bond>& bonds, const vector<Angle>& angles, const Xlink& xlink);
void readAtoms(const string& dataFileName);
void calculateBondLength(const vector <Bond>& bonds);
void testContactFile(const string& contactFileName);
void angleChecker(const string& dataFile);

int main(int argc, char* argv[])
{
	string dataFileName("4-MT-Mineralized-14.4576");
	string trajectoryFileName("dump.equilibrate");
	string outputDataFileName;
	string outputDirectory;
	string tempFileName("temp");
	string tclFileName("contact.tcl");
	string tclFileName2("contact2.tcl");
	string contactFileName("contacts.dat");
	string contactFileName2("contacts2.dat");
	string vmdCommand("vmd -dispdev text -e ");
	string vmdCommand2("vmd -dispdev text -e ");
	string healthStr;
	string cutoffStr, cutoffLoStr, percentStr;
	string logName("log.txt");
	ofstream log;
	vector<Bond> AGEbonds;
	vector<Bond> ENZbonds;
	vector<Bond> bonds;
	vector<Angle> AGEangles;
	vector<Angle> ENZangles;
	vector<Angle> angles;
	float cutoff = 16.70f;
	int initialBonds, numberOfCollagenAtoms, numberOfMolecules, beadsPerMolecule;
	int molsInRadius;
	float health;
	float radius, AGEradius;
	Xlink xlink = BTH;
	Coor centralAxis;
	char AGEAnswer;
	bool AGEInternal = false;
	bool debug = true;
	bool angleTest = false;
	bool mode1 = true;

	//testContactFile("contacts.dat");

	cout << "Mode 1?" << endl;
	char answer;
	cin >> answer;
	if (answer != 'y' && answer != 'Y') mode1 = false;

	if (argc == 1)
	{
		// Asking for inputs
		if (mode1) cout << "How much is the AGE percentage? (pmol/mgCollagen)" << endl;
		else cout << "How Healthy is this guy: (from 1 for diabetes to 4 for healthy)" << endl;
		cin >> health;
		cout << "Please enter the X and Y of the central axis:" << endl;
		cin >> centralAxis.x;
		cin >> centralAxis.y;
		cout << "Please enter the collagen microfibril radius and the AGE region radius:" << endl;
		cin >> radius;
		cin >> AGEradius;
		cout << "AGEs inside the AGE region:" << endl;
		cin >> AGEAnswer;
		if (AGEAnswer == 'y' || AGEAnswer == 'Y')
			AGEInternal = true;
	}
	else if (((argc == 2) ||(argc == 13)) && !mode1)
	{
		for (int i = 1; i < argc; ++i) {
			string arg = argv[i];
			if ((arg == "-h") || (arg == "--help")) {
				cout << "Usage: xLinker -i dataFileName -health healthNumber ";
				cout << "-i inputfileName -health healthValue -c centralAxisX centeralAxisY - r radius AGERadius - o outputDirectory " << endl;
				return 0;
			}
			else if (arg == "-i") {
				if (i + 1 < argc) { // Make sure we aren't at the end of argv!
					dataFileName = argv[++i];
				}
			}
			else if (arg == "-health") {
				if (i + 1 < argc) { // Make sure we aren't at the end of argv!
					health = atoi(argv[++i]);
				}
			}
			else if (arg == "-c") {
				if (i + 2 < argc) { // Make sure we aren't at the end of argv!
					centralAxis.x = atof(argv[++i]);
					centralAxis.x = atof(argv[++i]);
				}
			}
			else if (arg == "-r") {
				if (i + 2 < argc) { // Make sure we aren't at the end of argv!
					radius = atof(argv[++i]);
					AGEradius = atof(argv[++i]);
				}
			}

			else if (arg == "-o") {
				if (i + 1 < argc) { // Make sure we aren't at the end of argv!
					outputDirectory = argv[++i];
				}
			}
		}
	}

	stringstream ss;
	ss << health;
	ss >> healthStr;
	if (argc != 1)
	{
		outputDataFileName = outputDirectory + dataFileName + "-health-" + healthStr + ".dat";
		logName = outputDirectory + logName;
	}
	else
	{
		outputDataFileName = dataFileName + "-health-" + healthStr + ".dat";
	}
	ss.str("");
	//dataFileName += ".dat";
	// logger
	log.open(logName);
	if (mode1) log << "How much is the AGE percentage? (pmol/mgCollagen)" << endl;
	else log << "How Healthy is this guy: (from 1 for diabetes to 4 for healthy)" << endl;
	log << health << endl;
	log << "Please enter the X and Y of the central axis:" << endl;
	log << centralAxis.x << endl;
	log << centralAxis.y << endl;
	log << "Please enter the collagen microfibril radius and the AGE region radius:" << endl;
	log << radius << endl;
	log << AGEradius << endl;
	log << "AGEs inside the AGE region:" << endl;
	log << AGEAnswer << endl;




	beadsPerMolecule = measureBeadsPerMolecule(dataFileName, numberOfCollagenAtoms, numberOfMolecules);
	//writeTeminus(dataFileName, numberOfMolecules);
	//writeTeminus(trajectoryFileName, numberOfCollagenAtoms, beadsPerMolecule);
	molsInRadius = measureMolPerRadius(dataFileName, AGEradius, centralAxis);

	// Running VMD
	//writeAGETCLContact(tclFileName, trajectoryFileName, contactFileName, cutoff, 0.0f);
	//writeENZTCLContact(tclFileName2, dataFileName, contactFileName2, cutoff, 0.0f);
	vmdCommand  = vmdCommand  + contactFileName ;
	vmdCommand2 = vmdCommand2 + contactFileName2;
	//system(vmdCommand.c_str());
	//system(vmdCommand2.c_str());


	// X-Linking
	XLinker(health, xlink, dataFileName, contactFileName, contactFileName2, initialBonds, numberOfMolecules, molsInRadius, beadsPerMolecule, AGEbonds, ENZbonds,
		AGEangles, ENZangles, centralAxis, radius, AGEradius, log, mode1, AGEInternal);

	// Writing the data file
	if (xlink == BTH)
	{
		mergeBonds(bonds, AGEbonds, ENZbonds, debug);
		mergeAngles(angles, AGEangles, ENZangles);
	}
	writeDataFile(dataFileName, outputDataFileName, bonds, angles, xlink);

	if (angleTest) angleChecker(outputDataFileName);
	log.close();
	if (argc == 1)
	{
		cout << '\a';
		system("PAUSE");
	}
	else
	{ }
	return 0;
}

void writeAGETCLContact(const string& tclFileName, const string& trajectoryFileName,
	const string& contactFileName, const float& cutoff, const float& cutoffLo)
{
	ofstream tclFile(tclFileName);
	tclFile << "mol load lammpstrj " << trajectoryFileName << endl;
	tclFile << "set contactFile [open " << contactFileName << " w]" << endl;
	tclFile << "set type2 [atomselect top \"type 2\"]" << endl;
	tclFile << "set type3 [atomselect top \"type 3\"]" << endl;
	tclFile << "set contact2 [lindex [measure contacts " << cutoff << " $type2 $type3] 0]" << endl;
	tclFile << "set contact3 [lindex [measure contacts " << cutoff << " $type2 $type3] 1]" << endl;
	tclFile << "set contactSize [llength $contact2]" << endl;
	tclFile << "for {set i 0} {$i < $contactSize} {incr i} {" << endl;
	tclFile << "    set bondLength [measure bond [list [lindex $contact2 $i] [lindex $contact3 $i]]]" << endl;
	tclFile << "    if {$bondLength > " << cutoffLo << " } {" << endl;
	tclFile << "        puts $contactFile \"[lindex $contact2 $i] [lindex $contact3 $i]\"" << endl;
	tclFile << "    }" << endl;
	tclFile << "}" << endl;
	tclFile << "close $contactFile" << endl;
	tclFile << "exit" << endl;
}

void writeENZTCLContact(const string& tclFileName, const string& dataFileName,
	const string& contactFileName, const float& cutoff, const float& cutoffLo)
{
	string xyzFileName = dataFileName + "-ter.xyz";
	ofstream tclFile(tclFileName);
	tclFile << "mol load xyz " << xyzFileName << endl;
	tclFile << "set contactFile [open " << contactFileName << " w]" << endl;
	tclFile << "set type [atomselect top \"type 1\"]" << endl;\
	tclFile << "set contact1 [lindex [measure contacts " << cutoff << " $type] 0]" << endl;
	tclFile << "set contact2 [lindex [measure contacts " << cutoff << " $type] 1]" << endl;
	tclFile << "set contactSize [llength $contact1]" << endl;
	tclFile << "for {set i 0} {$i < $contactSize} {incr i} {" << endl;
	tclFile << "    set bondLength [measure bond [list [lindex $contact1 $i] [lindex $contact2 $i]]]" << endl;
	tclFile << "    if {$bondLength > " << cutoffLo << " } {" << endl;
	tclFile << "        puts $contactFile \"[lindex $contact1 $i] [lindex $contact2 $i]\"" << endl;
	tclFile << "    }" << endl;
	tclFile << "}" << endl;
	tclFile << "close $contactFile" << endl;
	tclFile << "exit" << endl;
}

void XLinker(const float& health, const Xlink& xlink,
	const string& dataFileName, const string& contactFileName,
	const string& contactFile2Name, const int&  initialBonds,
	int& numberOfMolecules, const int& molsInRadius,
	const int& beadsPerMolecule, vector<Bond>& AGEbonds,
	vector<Bond>& ENZbonds, vector<Angle>& AGEangles,
	vector<Angle>& ENZangles, const Coor& centralAxis,
	const float& radius, const float& AGEradius,
	ofstream& log, const bool& mode1, const bool& AGEInternal)
{
	vector<float> del;
	float healthTarget = health;
	int lastXNumber = 0;
	int XNumber = 0;
	int XTarget = 0;
	int counter = 0;
	int currentXTarget;
	switch(xlink)
	{
	case AGE:
		if (mode1) XTarget = (int)(health*numberOfMolecules / 10000.0f);
		else XTarget = (int)((6.0f - health)*numberOfMolecules / 10.0f);  //XTarget = (int)((6.0f-health)*molsInRadius/10.0f);
		break;
	case ENZ:
		if (mode1)
		{
			XTarget = 10000 + (int)(25.0f*(200.0f - health) / 24.0f);
			XTarget = 0; //(int)(XTarget*numberOfMolecules / 10000.0f);
		}
		else XTarget = (int)((1.0f + (1.0f / 3.0f)*(health - 1))*numberOfMolecules);
		break;
	case BTH:
		XLinker(health, AGE, dataFileName, contactFileName, contactFileName,
			initialBonds, numberOfMolecules, molsInRadius, beadsPerMolecule,
			AGEbonds, AGEbonds, AGEangles, AGEangles, centralAxis, radius,
			AGEradius, log, mode1, AGEInternal);
		XLinker(health, ENZ, dataFileName, contactFile2Name, contactFile2Name,
			initialBonds, numberOfMolecules, molsInRadius, beadsPerMolecule,
			ENZbonds, ENZbonds, ENZangles, ENZangles, centralAxis, radius,
			AGEradius, log, mode1, AGEInternal);
		return;
		break;
	}

	if (xlink == AGE)
		currentXTarget = 7 * XTarget /2 ;
	else
		currentXTarget = XTarget;

	while (XNumber != XTarget)
	{
		float dely = (float)(XNumber - XTarget);
		del.push_back(dely);
		counter++;
		if (abs(XTarget - XNumber) <= 5)
			currentXTarget += (XTarget - XNumber);
		else if (XTarget < XNumber)
			currentXTarget -= 5;
		else if (XTarget > XNumber)
			currentXTarget += 5;

		switch (xlink)
		{
		case AGE:
			XNumber = addAGE(dataFileName, contactFileName, currentXTarget, numberOfMolecules, beadsPerMolecule, AGEbonds, AGEangles, centralAxis, radius, AGEradius, AGEInternal);
			break;
		case ENZ:
			XNumber = addENZ(dataFileName, contactFile2Name, currentXTarget, numberOfMolecules, beadsPerMolecule,ENZbonds, ENZangles, xlink);
			break;
		}


		switch (xlink)
		{
		case AGE:
			cout << "Iteration " << counter << ": # of AGE is " << XNumber << " and currentXTarget is " << currentXTarget << " and XTarget is " << XTarget << endl;
			log << "Iteration " << counter << ": # of AGE is " << XNumber << " and currentXTarget is " << currentXTarget << " and XTarget is " << XTarget << endl;
			break;
		case ENZ:
			cout << "Iteration " << counter << ": # of ENZ is " << XNumber << " and currentXTarget is " << currentXTarget << " and XTarget is " << XTarget << endl;
			log << "Iteration " << counter << ": # of ENZ is " << XNumber << " and currentXTarget is " << currentXTarget << " and XTarget is " << XTarget << endl;
			break;
		}

		/*if (abs((XNumber - lastXNumber)) < 0.01f)
		{
			cout << "ABS" << endl;
			break;
		}*/
		if (counter == 70) break;  //70
		else if (abs(XNumber - XTarget) < 1 * XTarget / 100) break;
		else lastXNumber = XNumber;
	}
}

void writeTeminus(const string& dataFileName, int& numberOfMolecules)
{
	int numberOfCollagenAtoms;
	int numberOfCollagenMolecules;
	int beadsPerMolecule;
	string terminusXYZFileName = dataFileName + "-ter.xyz";
	string dataFileNameDat = dataFileName + ".dat";
	string line;
	ifstream dataFile(dataFileNameDat);
	ofstream terminusXYZFile(terminusXYZFileName);

	beadsPerMolecule = measureBeadsPerMolecule(dataFileName, numberOfCollagenAtoms, numberOfMolecules);
	numberOfCollagenMolecules = numberOfCollagenAtoms / beadsPerMolecule;

	terminusXYZFile << numberOfCollagenAtoms << endl;
	terminusXYZFile << "generated by MT" << endl;

	bool AtomsSection = false;
	int lineNumber = 0;
	while (!AtomsSection)
	{
		getline(dataFile, line);
		lineNumber++;
		int length = line.length();
		if (length >= 5)
			if (!strcmp(line.c_str(), "Atoms"))
			{
				getline(dataFile, line);
				AtomsSection = true;
			}
	}

	int atomID = 0;
	for (int i = 0; i < numberOfCollagenAtoms; i++)
	{
		getline(dataFile, line);
		atomID++;
		stringstream iss(line);
		int atomID, moleculeID, atomType;
		float x, y, z, d;
		vector<float> coors;
		iss >> atomID;
		iss >> moleculeID;
		iss >> atomType;
		iss >> d;
		iss >> x;
		iss >> y;
		iss >> z;
		coors.push_back(x);
		coors.push_back(y);
		coors.push_back(z);
		iss.clear();
		iss.str("");
		if ((atomID%beadsPerMolecule == 1) ||
			(atomID%beadsPerMolecule == 0))
			atomType = 1;
		else
			atomType = 2;
		XYZWriter(terminusXYZFileName, coors, atomType);
	}
}

void writeTeminus(const string& trajectoryFileName,
	const int& numberOfCollagenAtoms, const int& beadsPerMolecule)
{
	int numberOfCollagenMolecules;
	int numberOfBeads = 0;
	string terminusTrajFileName = trajectoryFileName + "-ter.lammpstrj";
	string trajFileNameDat = trajectoryFileName + ".lammpstrj";
	string line;
	ifstream trajFile(trajFileNameDat);
	ofstream terminusTrajFile(terminusTrajFileName);

	numberOfCollagenMolecules = numberOfCollagenAtoms / beadsPerMolecule;

	int lineNumber = 0;
	int id, type;
	float x, y, z;
	while (getline(trajFile, line))
	{
		terminusTrajFile << line << endl;
		lineNumber++;
		int lineLength = line.size();

		if (lineLength == 21)
		{
			if (!strcmp(line.c_str(), "ITEM: NUMBER OF ATOMS"))
			{
				getline(trajFile, line);
				terminusTrajFile << line << endl;
				stringstream iss(line);
				iss >> numberOfBeads;
				iss.clear();
				iss.str("");
			}
		}
		else if (lineLength == 28)
		{
			if (!strcmp(line.c_str(), "ITEM: ATOMS id type xs ys zs"))
			{
				for (int i = 0; i < numberOfBeads; i++)
				{
					getline(trajFile, line);
					stringstream iss(line);
					iss >> id;
					iss >> type;
					iss >> x;
					iss >> y;
					iss >> z;

					if ((id%beadsPerMolecule == 1) ||
						(id%beadsPerMolecule == 0))
						type = 1;
					else
						type = 2;
					if (id >= numberOfCollagenAtoms)
						type = 3;
					terminusTrajFile << id << " "
						<< type << " "
						<< x << " "
						<< y << " "
						<< z << endl;
					iss.clear();
					iss.str("");
				}
			}
		}
	}
}

void XYZWriter(const string& xyzFileName, const vector<float> coors, const int& atomType)
{
	float x, y, z;
	ofstream file(xyzFileName, ios_base::app);
	x = coors[0];
	y = coors[1];
	z = coors[2];

	file << " ";
	file.width(8);
	file << fixed << left << atomType;
	file.width(10);
	file << fixed << setprecision(6) << right << x;
	file << "      ";
	file.width(10);
	file << fixed << setprecision(6) << right << y;
	file << "      ";
	file.width(10);
	file << fixed << setprecision(6) << right << z;
	file << endl;
}

int addAGE(const string& dataFileName, const string& contactFileName,
	const int& XTarget, int& numberOfMolecules, const int&  beadsPerMolecule,
	vector<Bond>& bonds, vector<Angle>& angles, const Coor& centralAxis,
	const float& radius, const float& AGEradius, const bool& AGEInternal)
{
	bonds.clear();
	angles.clear();
	string dataFileNameDat = dataFileName + ".dat";
	ifstream contactFile(contactFileName);
	ifstream dataFile(dataFileNameDat);
	int numberOfCollagenAtoms;
	string line;
	int numberOfCrosslinks = 0;
	int realNumberOfCrosslinks = 0;
	float zLow, zHi;
	int numberOfZSections = 0;
	int AGEinSection = 10;
	vector<Bond> contacts;
	vector<int> atomSections;
	vector<float> zSectionLow;
	vector<bool> atomFlags;
	vector<int> zSectionFlags;
	vector<bool> moleculeFlags;
	vector<Coor> coors;
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(0, 100000);
	vector<double> randomNumbers;
	int percentOfChosenCrosslinks;
	int angleTypes;

	if (XTarget < 10) AGEinSection = 1;

	measureBeadsPerMolecule(dataFileName, numberOfCollagenAtoms, numberOfMolecules);

	for (int i = 0; i < numberOfMolecules; i++)
		moleculeFlags.push_back(false);

	random_device random_device;
	mt19937 random_engine(random_device());
	std::uniform_int_distribution<int> distribution(1, numberOfMolecules);

	for (int i = 0; i < XTarget; i++)
	{
		int randomMolecule = distribution(random_engine);
		moleculeFlags[randomMolecule-1] = true;
	}

	while (getline(contactFile, line))
	{
		stringstream iss(line);
		Bond contact;
		int one, two;
		int type = 2;
		iss >> one;
		iss >> two;
		contact.first = one;
		contact.second = two;
		contact.type = type;
		contacts.push_back(contact);
		numberOfCrosslinks++;
	}

	for (int i = 0; i < numberOfCrosslinks; i++)
		randomNumbers.push_back(dis(gen));
	percentOfChosenCrosslinks = (100000*XTarget) / numberOfCrosslinks;

	int lineNumber = 0;
	int flag = true;
	while (flag)
	{
		getline(dataFile, line);
		lineNumber++;
		int lineSize = line.length();

		if (lineNumber == 15)
		{
			numberOfZSections = XTarget/AGEinSection;
			stringstream iss("");
			iss << line;
			iss >> zLow;
			iss >> zHi;
			for (int i = 0; i < numberOfZSections; i++)
			{
				float sectionLow = zLow + i*(zHi - zLow) / numberOfZSections;
				zSectionFlags.push_back(0);
				zSectionLow.push_back(sectionLow);
			}
		}
		else if (lineNumber == 11)
		{
			stringstream iss("");
			iss << line;
			iss >> angleTypes;
			iss.clear();
			iss.str("");
		}
		else if (lineSize >= 5)
		{
			if (!strcmp(line.c_str(), "Atoms"))
			{
				getline(dataFile, line); // Skipping the empty line
				flag = false;
				for (int i = 0; i < numberOfCollagenAtoms; i++)
				{
					int dumInt;
					float dumFloat;
					float xVal,yVal,zVal;
					Coor coorVal;
					int sectionNumber;
					getline(dataFile, line);
					stringstream iss("");
					iss << line;
					iss >> dumInt;
					iss >> dumInt;
					iss >> dumInt;
					iss >> dumFloat;
					iss >> xVal;
					iss >> yVal;
					iss >> zVal;
					coorVal.x = xVal;
					coorVal.y = yVal;
					coors.push_back(coorVal);
					sectionNumber = (int)((zVal - zLow) / ((zHi - zLow) / numberOfZSections));
					atomSections.push_back(sectionNumber);
				}
			}
		}

	}

	for (int i = 0; i < numberOfCrosslinks; i++)
	{
		int firstMoleculeID = contacts[i].first / beadsPerMolecule;
		int secondMoleculeID = contacts[i].second / beadsPerMolecule;
		int firstSection  = atomSections[contacts[i].first - 1];
		int secondSection = atomSections[contacts[i].second - 1];
		/*if (moleculeFlags[firstMoleculeID] &&
			moleculeFlags[secondMoleculeID]) &&*/
	    /*if ((zSectionFlags[firstSection]  < AGEinSection) &&
			(zSectionFlags[secondSection] < AGEinSection))*/
		if (randomNumbers[i] < percentOfChosenCrosslinks)
		{
			float xValFirst, yValFirst;
			float xValSecond, yValSecond;
			float firstR, secondR;
			float axisX, axisY;
			xValFirst = coors[contacts[i].first - 1].x;
			yValFirst = coors[contacts[i].first - 1].y;
			xValSecond = coors[contacts[i].second - 1].x;
			yValSecond = coors[contacts[i].second - 1].y;
			axisX = centralAxis.x;
			axisY = centralAxis.y;
			firstR  = sqrt((axisX -  xValFirst)*(axisX -  xValFirst) + (axisY -  yValFirst)*(axisY -  yValFirst));
			secondR = sqrt((axisX - xValSecond)*(axisX - xValSecond) + (axisY - yValSecond)*(axisY - yValSecond));
			if ( (  (firstR > AGEradius)  &&  (secondR > AGEradius)  && !AGEInternal) ||
				 (  (firstR < AGEradius)  &&  (secondR < AGEradius)  &&  AGEInternal)  )
			{
				realNumberOfCrosslinks++;
				bonds.push_back(contacts[i]);
				moleculeFlags[firstMoleculeID] = false;
				moleculeFlags[secondMoleculeID] = false;
				zSectionFlags[atomSections[contacts[i].first - 1]]++;
				zSectionFlags[atomSections[contacts[i].second - 1]]++;

				Angle myAngle;
				myAngle.type = angleTypes + 1;
				if ((contacts[i].first%beadsPerMolecule) &&
					(contacts[i].first - 1 != contacts[i].second))
				{
					myAngle.first = contacts[i].first - 1;
					myAngle.second = contacts[i].first;
					myAngle.third = contacts[i].second;
					angles.push_back(myAngle);
				}
				if ((contacts[i].first%beadsPerMolecule != beadsPerMolecule - 1) &&
					(contacts[i].first + 1 != contacts[i].second))
				{
					myAngle.first = contacts[i].first + 1;
					myAngle.second = contacts[i].first;
					myAngle.third = contacts[i].second;
					angles.push_back(myAngle);
				}
				if ((contacts[i].second%beadsPerMolecule) &&
					(contacts[i].second - 1 != contacts[i].first))
				{
					myAngle.first = contacts[i].second - 1;
					myAngle.second = contacts[i].second;
					myAngle.third = contacts[i].first;
					angles.push_back(myAngle);
				}
				if ((contacts[i].second%beadsPerMolecule != beadsPerMolecule - 1) &&
					(contacts[i].first != contacts[i].second + 1))
				{
					myAngle.first = contacts[i].second + 1;
					myAngle.second = contacts[i].second;
					myAngle.third = contacts[i].first;
					angles.push_back(myAngle);
				}
			}
		}
	}


	return realNumberOfCrosslinks;
}

int addENZ(const string& dataFileName, const string& contactFileName,
	const int& XTarget, int& numberOfMolecules, const int&  beadsPerMolecule,
	vector<Bond>& bonds, vector<Angle>& angles,
	const Xlink& xlink)
{
	bonds.clear();
	angles.clear();
	string dataFileNameDat = dataFileName + ".dat";
	ifstream contactFile(contactFileName);
	ifstream dataFile(dataFileNameDat);
	string line;
	int numberOfCrosslinks = 0;
	int realNumberOfCrosslinks = 0;
	vector<Bond> contacts;
	vector<int> atomSections;
	vector<float> zSectionLow;
	vector<bool> atomFlags;
	vector<bool> zSectionFlags;
	vector<bool> moleculeFlags;
	vector<Coor> coors;
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(0, 100000);
	vector<double> randomNumbers;
	int angleTypes = 0;
	int percentOfChosenCrosslinks;

	// Checking if the files are found:
	if (contactFile.fail() || dataFile.fail())
	{
		cout << "Either " << contactFileName << " or " << dataFileName << " is not found" << endl;
		return 0;
	}

	for (int i = 0; i < numberOfMolecules; i++)
		moleculeFlags.push_back(false);

	random_device random_device;
	mt19937 random_engine(random_device());
	std::uniform_int_distribution<int> distribution(1, numberOfMolecules);

	for (int i = 0; i < XTarget; i++)
	{
		int randomMolecule = distribution(random_engine);
		moleculeFlags[randomMolecule - 1] = true;
	}

	while (std::getline(contactFile, line))
	{
		stringstream iss(line);
		Bond contact;
		int one, two;
		int type;
		switch (xlink)
		{
		case ENZ:
			type = 2;
			break;
		case BTH:
			type = 3;
			break;
		}
		iss >> one;
		iss >> two;
		contact.first = one;
		contact.second = two;
		contact.type = type;
		contacts.push_back(contact);
		numberOfCrosslinks++;
	}

	int lineNumber = 0;
	bool flag = true;
	while (flag)
	{
		getline(dataFile, line);
		lineNumber++;

		if (lineNumber == 11)
		{
			stringstream iss("");
			iss << line;
			iss >> angleTypes;
			iss.clear();
			iss.str("");
			flag = false;
		}
	}

	for (int i = 0; i < numberOfCrosslinks; i++)
		randomNumbers.push_back(dis(gen));
	percentOfChosenCrosslinks = (100000 * XTarget) / numberOfCrosslinks;

	for (int i = 0; i < numberOfCrosslinks; i++)
	{
		int firstMoleculeID = contacts[i].first / beadsPerMolecule;
		int secondMoleculeID = contacts[i].second / beadsPerMolecule;
		if (randomNumbers[i] < percentOfChosenCrosslinks)
		{
			realNumberOfCrosslinks++;
			bonds.push_back(contacts[i]);
			realNumberOfCrosslinks;
			Angle myAngle;
			myAngle.type = angleTypes + 1;
			if ((contacts[i].first%beadsPerMolecule) &&
				(contacts[i].first - 1 != contacts[i].second))
			{
				myAngle.first = contacts[i].first - 1;
				myAngle.second = contacts[i].first;
				myAngle.third = contacts[i].second;
				angles.push_back(myAngle);
			}
			if ((contacts[i].first%beadsPerMolecule != beadsPerMolecule - 1) &&
				(contacts[i].first + 1 != contacts[i].second))
			{
				myAngle.first = contacts[i].first + 1;
				myAngle.second = contacts[i].first;
				myAngle.third = contacts[i].second;
				angles.push_back(myAngle);
			}
			if ((contacts[i].second%beadsPerMolecule) &&
				(contacts[i].second - 1 != contacts[i].first))
			{
				myAngle.first = contacts[i].second - 1;
				myAngle.second = contacts[i].second;
				myAngle.third = contacts[i].first;
				angles.push_back(myAngle);
			}
			if ((contacts[i].second%beadsPerMolecule != beadsPerMolecule - 1) &&
				(contacts[i].second + 1 != contacts[i].first))
			{
				myAngle.first = contacts[i].second + 1;
				myAngle.second = contacts[i].second;
				myAngle.third = contacts[i].first;
				angles.push_back(myAngle);
			}
		}
	}
	return realNumberOfCrosslinks;
}

int measureBeadsPerMolecule(const string& dataFileName, int& numberOfCollagenAtoms, int& numberOfMolecules)
{
	int beadsPerMolecule = 0;
	int lineNumber = 0;
	int numberOfAtomTypes = 0;
	bool atomsSection = false;
	string line;
	string cDataFileName = dataFileName + ".dat";
	ifstream dataFile(cDataFileName);

	// Checking if the file exist
	if (dataFile.fail())
	{
		cout << cDataFileName << " is not found" << endl;
		return 0;
	}

	numberOfCollagenAtoms = 0;
	bool flag = true;
	while (flag)
	{
		getline(dataFile, line);
		lineNumber++;
		if (lineNumber == 9)
		{
			stringstream iss(line);
			iss >> numberOfAtomTypes;
			iss.clear();
			iss.str("");
		}

		if (atomsSection)
		{
			stringstream iss(line);
			int atomID, moleculeID, atomType;
			iss >> atomID;
			iss >> moleculeID;
			iss >> atomType;
			iss.str("");
			numberOfCollagenAtoms++;
			if (moleculeID == 1)
			{
				beadsPerMolecule++;
			}
			else if (atomType == numberOfAtomTypes)
			{
				numberOfCollagenAtoms--;
				flag = false;
				break;
			}
			else
			{
				numberOfMolecules = moleculeID;
			}
		}
		else if (lineNumber == 9)
		{
			stringstream iss3(line);
			iss3 >> numberOfAtomTypes;
		}
		else if (line.length() >= 5)
			if (!strcmp(line.c_str(), "Atoms"))
			{
				getline(dataFile, line);
				atomsSection = true;
			}
	}
	dataFile.close();
	return beadsPerMolecule;
}

int measureMolPerRadius(const string& dataFileName, const float& AGEradius, const Coor& centralAxis)
{
	int beadsPerMolecule;
	int numberOfCollagenAtoms;
	int numberOfMolecules;
	int lineNumber = 0;
	string line;
	string cDataFileName = dataFileName + ".dat";
	ifstream dataFile(cDataFileName);
	int beadsInRadius = 0;
	int molsInRadius;

	beadsPerMolecule = measureBeadsPerMolecule(dataFileName, numberOfCollagenAtoms, numberOfMolecules);



	bool atomsSection = false;
	while (!atomsSection)
	{
		getline(dataFile, line);
		int lineLength = line.size();

		if (line.length() >= 5)
			if (!strcmp(line.c_str(), "Atoms"))
			{
				getline(dataFile, line);
				atomsSection = true;
			}
	}

	for (int i = 0; i < numberOfCollagenAtoms; i++)
	{
		getline(dataFile, line);
		int id, molId, type;
		float q, x, y, z, R;
		stringstream iss(line);
		iss >> id;
		iss >> molId;
		iss >> type;
		iss >> q;
		iss >> x;
		iss >> y;
		iss >> z;
		iss.clear();
		iss.str("");

		R = sqrt((x - centralAxis.x)*(x - centralAxis.x) + (y - centralAxis.y)*(y - centralAxis.y));
		if (R > AGEradius)
			beadsInRadius++;
	}

	molsInRadius = beadsInRadius / beadsPerMolecule;

	return molsInRadius;
}

void mergeBonds(vector<Bond>& bonds, const vector<Bond>& AGEbonds, const vector<Bond>& ENZbonds, const bool debug)
{
	int AGELength = AGEbonds.size();
	int ENZLength = ENZbonds.size();
	vector<Bond> ENZbondsCopy = ENZbonds;

	for (int i = 0; i < AGELength; i++)
		bonds.push_back(AGEbonds[i]);
	for (int i = 0; i < ENZLength; i++)
	{
		ENZbondsCopy[i].type = 3;
		bonds.push_back(ENZbondsCopy[i]);
	}

	if (debug)
	{
		ofstream bondLength("BondLength.dat");
		readAtoms("4-MT-Mineralized-9.05211");
		int numberOfBonds = bonds.size();
		for (int i = 0; i < numberOfBonds; i++)
		{
			int id1, id2, type;
			float x1, y1, z1, x2, y2, z2;
			float l;
			id1 = bonds[i].first;
			id2 = bonds[i].second;
			type = bonds[i].type;
			x1 = atoms[id1].x;
			x2 = atoms[id2].x;
			y1 = atoms[id1].y;
			y2 = atoms[id2].y;
			z1 = atoms[id1].z;
			z2 = atoms[id2].z;
			l = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
			bondLength << id1 << "," <<
				id2 << "," << type << " " <<
				l << endl;
		}
	}
}

void mergeAngles(vector<Angle>& angles, const vector<Angle>& AGEangles, const vector<Angle>& ENZangles)
{
	int AGELength = AGEangles.size();
	int ENZLength = ENZangles.size();
	vector<Angle> ENZanglesCopy = ENZangles;

	for (int i = 0; i < AGELength; i++)
		angles.push_back(AGEangles[i]);
	for (int i = 0; i < ENZLength; i++)
	{
		ENZanglesCopy[i].type += 1;
		angles.push_back(ENZanglesCopy[i]);
	}
}

void writeDataFile(const string& dataFileName, const string& outputDataFileName,
	const vector<Bond>& bonds, const vector<Angle>& angles, const Xlink& xlink)
{
	string dataFileNameDat = dataFileName + ".dat";
	ifstream dataFile(dataFileNameDat);
	ofstream outputDataFile(outputDataFileName);
	string line;
	int oldNumberOfBonds, oldNumberOfAngles;
	int numberOfNewBonds = bonds.size();
	int numberOfNewAngles = angles.size();
	int angleTypes;

	if (dataFile.fail())
	{
		cout << dataFileName << " is not found" << endl;
		return;
	}


	int lineNumber = 0;
	bool flag = true;
	while (flag)
	{
		getline(dataFile, line);
		lineNumber++;
		int lineLength = line.length();
		if (lineNumber == 4)
		{
			int bonds;
			stringstream iss("");
			iss << line;
			iss >> oldNumberOfBonds;
			bonds = numberOfNewBonds + oldNumberOfBonds;
			outputDataFile << bonds << " bonds" << endl;
		}
		else if (lineNumber == 5)
		{
			int angles;
			stringstream iss("");
			iss << line;
			iss >> oldNumberOfAngles;
			angles = numberOfNewAngles + oldNumberOfAngles;
			outputDataFile << angles << " angles" << endl;
		}
		else if (lineNumber == 10)
		{
			switch (xlink)
			{
			case AGE:
			case ENZ:
				outputDataFile << "2 bond types" << endl;
				break;
			case BTH:
				outputDataFile << "3 bond types" << endl;
				break;
			default:
				break;
			}
		}
		else if (lineNumber == 11)
		{
			stringstream iss("");
			iss << line;
			iss >> angleTypes;
			iss.clear();
			iss.str("");
			switch (xlink)
			{
			case AGE:
			case ENZ:
				outputDataFile << angleTypes + 1;
				outputDataFile << " angle types" << endl;
				break;
			case BTH:
				outputDataFile << angleTypes + 2;
				outputDataFile << " angle types" << endl;
				break;
			default:
				break;
			}

		}
		else if (lineLength >= 11)
		{
			outputDataFile << line << endl;
			if (!strcmp(line.c_str(), "Bond Coeffs"))
			{
				getline(dataFile, line);
				outputDataFile << line << endl;  // empty line
				getline(dataFile, line);
				outputDataFile << "  1 8.565 14.00 62.780 18.20 00.000 21.00 21.00" << endl;
				switch (xlink)
				{
				case AGE:
					outputDataFile << "  2 0.440 10.27 10.605 14.04 46.790 16.29 17.63" << endl;
					break;
				case ENZ:
					outputDataFile << "  2 2.140 11.62 11.626 12.50 44.260 12.74 13.39" << endl;
					break;
				case BTH:
					outputDataFile << "  2 0.440 10.27 10.605 14.04 46.790 16.29 17.63" << endl;
					outputDataFile << "  3 2.140 11.62 11.626 12.50 44.260 12.74 13.39" << endl;
					break;
				}
			}
			else if (!strcmp(line.c_str(), "Angle Coeffs"))
			{
				getline(dataFile, line);
				outputDataFile << line << endl;
				for (int i = 0; i < angleTypes; i++)
				{
					int id, theta;
					float k;
					getline(dataFile, line);
					stringstream iss(line);
					iss >> id;
					iss >> k;
					iss >> theta;
					outputDataFile << "  " << i + 1 << " " << k << " " << theta << endl;
					iss.clear();
					iss.str("");
				}

				switch (xlink)
				{
				case AGE:
					outputDataFile << "  " << angleTypes + 1 <<
						" 22 90 # The actual K value is much higher; use the shake constraint instead" << endl;
					break;
				case ENZ:
					outputDataFile << "  " << angleTypes + 1 <<
						" 22 90 # The actual K value is much higher; use the shake constraint instead" << endl;
					break;
				case BTH:
					outputDataFile << "  " << angleTypes + 1 <<
						" 22 90 # The actual K value is much higher; use the shake constraint instead" << endl;
					outputDataFile << "  " << angleTypes + 2 <<
						" 22 90 # The actual K value is much higher; use the shake constraint instead" << endl;
					break;
				}
			}
		}
		else if (lineLength >= 5)
		{
			outputDataFile << line << endl;
			if (!strcmp(line.c_str(), "Bonds"))
			{
				getline(dataFile, line); // empty line
				outputDataFile << line << endl;
				flag = false;
				for (int i = 0; i < oldNumberOfBonds; i++)
				{
					getline(dataFile, line);
					outputDataFile << line << endl;
				}
			}
			else if (!strcmp(line.c_str(), "Atoms"))
			{}
		}
		else
			outputDataFile << endl;
	}

	for (int i = 0; i < numberOfNewBonds; i++)
		outputDataFile << "  " << oldNumberOfBonds + i + 1
			<< " " << bonds[i].type
			<< " " << bonds[i].first + 1
			<< " " << bonds[i].second + 1
			<< endl; // VMD writes indexes which are zero-based, while LAMMPS uses one-based serials.


	for (int i = 0; i < 3 + oldNumberOfAngles; i++)
	{
		getline(dataFile, line);
		outputDataFile << line << endl;
	}

	for (int i = 0; i < numberOfNewAngles; i++)
	{
		int first, second, third, type;
		type = angles[i].type;
		first = angles[i].first;
		second = angles[i].second;
		third = angles[i].third;
		outputDataFile << "  " << oldNumberOfAngles + i + 1
			<< " " << type
			<< " " << first + 1
			<< " " << second + 1
			<< " " << third + 1
			<< endl; // VMD writes indexes which are zero-based, while LAMMPS uses one-based serials.
	}

	while (getline(dataFile, line))
		outputDataFile << line << endl;
}

void readAtoms(const string& dataFileName)
{
	int numberOfCollagenAtoms;
	int numberOfMolecules;
	string dataFileNameDat = dataFileName + ".dat";
	string line;
	ifstream dataFile(dataFileNameDat);
	measureBeadsPerMolecule(dataFileName, numberOfCollagenAtoms, numberOfMolecules);

	bool AtomsSection = false;
	int lineNumber = 0;
	while (!AtomsSection)
	{
		getline(dataFile, line);
		lineNumber++;
		int length = line.length();
		if (length >= 5)
			if (!strcmp(line.c_str(), "Atoms"))
			{
				getline(dataFile, line);
				AtomsSection = true;
			}
	}

	int atomID = 0;
	for (int i = 0; i < numberOfCollagenAtoms; i++)
	{
		getline(dataFile, line);
		atomID++;
		stringstream iss(line);
		int atomID, moleculeID, atomType;
		float x, y, z, d;
		Atom atom;
		iss >> atomID;
		iss >> moleculeID;
		iss >> atomType;
		iss >> d;
		iss >> x;
		iss >> y;
		iss >> z;
		atom.id = atomID;
		atom.x = x;
		atom.y = y;
		atom.z = z;
		atoms.push_back(atom);
		iss.clear();
		iss.str("");
	}
	dataFile.close();
}

void calculateBondLength(const vector <Bond>& bonds)
{
	int firstID, secondID;
	float x1, y1, z1, x2, y2, z2;
	float distance;
	int numberOfBonds = bonds.size();
	ofstream bondFile("bonds.txt");
	for (int i = 0; i < numberOfBonds; i++)
	{
		firstID = bonds[i].first;
		secondID = bonds[i].second;
		x1 = atoms[firstID].x;
		y1 = atoms[firstID].y;
		z1 = atoms[firstID].z;
		x2 = atoms[secondID].x;
		y2 = atoms[secondID].y;
		z2 = atoms[secondID].z;
		distance = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
		bondFile << firstID << "," <<
			       secondID << " " <<
			       distance << endl;
	}
}

void testContactFile(const string& contactFileName)
{
	readAtoms("4-MT-Mineralized-9.05211");
	ifstream contactFile(contactFileName);
	vector<Bond> bonds;
	int first, second,type;
	string line;
	int lineNumber = 0;
	while (getline(contactFile, line))
	{
		lineNumber++;
		if (lineNumber == 2)
			cout << line << endl;
		if (lineNumber == 4)
			cout << line << endl;
		Bond bond;
		stringstream iss(line);
		iss >> first;
		iss >> second;
		iss.clear();
		iss.str("");
		type = 1;
		bond.first = first;
		bond.second = second;
		bond.type = type;
		bonds.push_back(bond);
	}
	calculateBondLength(bonds);
}

void angleChecker(const string& dataFileName)
{
	ifstream dataFile(dataFileName);
	string line;
	vector<Bond> bonds;
	vector<Angle> angles;
	int lineNumber = 0;
	int numBonds;
	int numAngles;
	stringstream iss;

	while (getline(dataFile, line))
	{
		int id1, id2, id3;
		int type, serial;
		lineNumber++;
		if (lineNumber == 4)
		{
			iss << line;
			iss >> numBonds;
			iss.clear();
			iss.str("");
		}
		else if (lineNumber == 5)
		{
			iss << line;
			iss >> numAngles;
			iss.clear();
			iss.str("");
		}
		else if (!strcmp(line.c_str(), "Bonds"))
		{
				getline(dataFile, line);
				for (int i = 0; i < numBonds; i++)
				{
					getline(dataFile, line);
					iss << line;
					iss >> serial;
					iss >> type;
					iss >> id1;
					iss >> id2;
					iss.clear();
					iss.str("");
					Bond myBond;
					myBond.first = id1;
					myBond.second = id2;
					myBond.type = type;
					bonds.push_back(myBond);
				}
				getline(dataFile, line);
				getline(dataFile, line);
				getline(dataFile, line);
				for (int i = 0; i < numAngles; i++)
				{
					getline(dataFile, line);
					iss << line;
					iss >> serial;
					iss >> type;
					iss >> id1;
					iss >> id2;
					iss >> id3;
					iss.clear();
					iss.str("");
					Angle myAngle;
					myAngle.first = id1;
					myAngle.second = id2;
					myAngle.third = id3;
					myAngle.type = type;
					angles.push_back(myAngle);
				}
		}
	}

	int counter = 0;
	for ( ; counter < numAngles; counter++)
	{
		cout << "Checking angle " << counter << endl;
		int j = 0;
		int k = 0;
		int id1 = angles[counter].first;
		int id2 = angles[counter].second;
		int id3 = angles[counter].third;
		int type = angles[counter].type;
		if (type == 1) continue;
		for ( ; j < numBonds; j++)
			if ((bonds[j].first == id1 && bonds[j].second == id2)  ||
				(bonds[j].first == id2 && bonds[j].second == id1))
				break;
		for (; k < numBonds; k++)
			if ((bonds[k].first == id2 && bonds[k].second == id3) ||
				(bonds[k].first == id3 && bonds[k].second == id2))
				break;
		if (j == numBonds)
		{
			cout << "Bond " << id1 << "," << id2 << " not found!!" << endl;
			break;
		}
		if (k == numBonds)
		{
			cout << "Bond " << id2 << "," << id3 << " not found!!" << endl;
			break;
		}
	}

	if (counter == numAngles)
		cout << "Angles passed the test !!!" << endl;
	else
		cout << "Angles did not pass the test !!!" << endl;
}
