class ContainerConfigurer {
    static String process(def params) {
        def container = params?.rn_container?.toString()?.trim()

        if (!container) {
            return container
        }

        if (!container.startsWith('file:///')) {
            throw new IllegalArgumentException(
                "SC_BLIPPER_CONTAINER_PATH must use an absolute file:/// URI: ${container}"
            )
        }

        def image = new File(new URI(container))
        if (!image.isFile() || !image.canRead()) {
            throw new IllegalArgumentException(
                "SC_BLIPPER_CONTAINER_PATH is not a readable file on the Nextflow driver: ${image}"
            )
        }

        container
    }
}
