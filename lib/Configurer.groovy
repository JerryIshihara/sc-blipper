class Configurer {
    static Map process(def params) {
        [
            cnmf: CnmfConfigurer.process(params),
            scvi: ScviConfigurer.process(params)
        ]
    }
}
