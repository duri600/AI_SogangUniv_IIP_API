class ProcBuffers
{
private:
	
public:
	ProcBuffers();
	~ProcBuffers();
	static void Process(double** input, int type);
	void Training(int type);
};