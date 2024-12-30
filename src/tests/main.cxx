#define NOMINMAX

#include <cstdlib>

#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>

class testRunListener : public Catch::EventListenerBase {
public:
    using Catch::EventListenerBase::EventListenerBase;

    void testRunStarting(Catch::TestRunInfo const&) override {
        
//        auto str_env = [](const char* var) -> std::string {
//            const char* chr = std::getenv(var);
//            if (chr == nullptr){ return ""; }
//            return chr;
//        };
//        std::cout << "RPPREFIX: " << str_env("RPPREFIX") << std::endl;
//        std::cout << "RESOURCES: " << str_env("RESOURCES") << std::endl;
//        if (!std::filesystem::exists(str_env("RPPREFIX"))){
//            std::cerr << "RPPREFIX variable was not set to a valid path; terminating" << std::endl;
//            std::exit(100);
//        }
//        if (!std::filesystem::exists(str_env("RESOURCES"))){
//            std::cerr << "RESOURCES variable was not set to a valid path" << std::endl;
//        }
//        
//        auto load_method = AbstractSharedLibraryWrapper::load_method::LOAD_LIBRARY;
//        
//        namespace fs = std::filesystem;
//        fs::path target(fs::path(std::getenv("RPPREFIX")) / fs::path("librefprop.dylib"));
//        std::unique_ptr<NativeSharedLibraryWrapper> RP;
//        RP.reset(new NativeSharedLibraryWrapper(target.string(), load_method));
    }
};
CATCH_REGISTER_LISTENER(testRunListener)
